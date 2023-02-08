
#summarises somatic SNVs and CNV calls into subclone evolution over samples
#and determines which subclones are subclones of which other subclones.

#' Combines variants and CNAs to clonal evolution
#'
#' @param variants variants: The variants.
#' @param cnvs cnvs: The copy number calls.
#' @param timeSeries names list of vectors: samples to be analysed together, named by individual.
#' @param normals named boolean vector: which sameples (names of vector) are normal.
#' @param genome character: the genome.
#' @param Rdirectory character: The save directory.
#' @param plotDirectory character: the directory to plot to.
#' @param cpus integer: maximum number of parallel processes. Default 1.
#' @param forceRedo boolen: if redoing calculations even if saved data is available. Default FALSE.
#'
#' @export
#'
#' @details This function calls VEP on the output from outputSomaticVariants. For this, VEP needs to be callable by system('vep').
getStories = function(variants, cnvs, timeSeries, normals, genome, cloneDistanceCut=-qnorm(0.01),
  Rdirectory, plotDirectory, cpus=1, forceRedo=F, manualStoryMerge=F, correctReferenceBias=T, rareGermline=T, maxStories = 3000) {
  stories = list()
  saveFile = paste0(Rdirectory, '/stories.Rdata')
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading stories.\n')
    load(file=saveFile)
    return(stories)
  }

  #add theoretical errors to the CNA calls
  #this accounts for things such as 100% AAB = 50% AAAB.
  #removed (for now), as the consistency and dodgyness fixes handles the same problem
  #and fixing it twice makes it too conservative.
  #cnvs = addTheoreticalCNVerrors(cnvs)
  
  if ( length(timeSeries) > 0 ) {
    for ( i in 1:length(timeSeries) ) {
      ts = timeSeries[[i]]
      name = names(timeSeries)[i]
      qs = variants$variants[ts]
      catLog('\nTracking clonal evolution in ', name, '..', sep='')
      
      #select somatic SNPs
      if ( any(!normals[ts]) ) {
        somaticMx = do.call(cbind, lapply(qs[!normals[ts]], function(q) q$somaticP > 0.95))
        somatic = apply(somaticMx, 1, any)
        catLog('Found ', sum(somatic), ' high quality somatic SNVs to track.\n', sep='')
      }
      else
        somatic = rep(FALSE, nrow(qs[[1]]))

      somaticQs = lapply(qs, function(q) q[somatic,])

      #switch to effective coverage, to not overestimate the accuracy of high-coverage SNPs
      #not too low though, keep at least 100.
      mC = max(100, getMaxCov())
      somaticQs = lapply(somaticQs, function(q) {
        d = q$cov
        effectiveCov = round(d*(1 + d/mC)/(1 + d/mC + d^2/mC^2))
        effectiveVar = round(q$var/q$cov*effectiveCov)
        q$var = effectiveVar
        q$cov = effectiveCov
        q$ref = q$cov - q$var
        return(q)
      })
      
      #set clonality of SNPs from frequency and local CNV, return stories
      catLog('SNVs..\n')
      snpStories = findSNPstories(somaticQs, cnvs[ts], normals[ts], filter=T)
      if ( nrow(snpStories) > maxStories ) {
        catLog('restricting to ', maxStories, ' random somatic SNVs in the interest of time. Rest will be added back in after clustering.\n', sep='')
        snpStories = snpStories[sample(1:nrow(snpStories), maxStories),]
      }
      anchorSNVs = snpStories
      catLog('Keeping', nrow(snpStories), 'SNV stories.\n')
      
      #combine CNV calls over sample into stories
      catLog('Tracking clonal evolution in CNVs..')
      cnvStories = getCNVstories(cnvs[ts], normals[ts], genome, filter=T)
      anchorCNAs = cnvStories
      catLog('Keeping', nrow(cnvStories), 'CNV stories.\n')

      #merge SNV and CNA stories.
      catLog('merge events into clones..')
      allStories = data.frame(row.names=c('germline', rownames(snpStories), rownames(cnvStories)), stringsAsFactors=F)
      allStories$x1 = c(NA, snpStories$x1, cnvStories$x1)
      allStories$x2 = c(NA, snpStories$x2, cnvStories$x2)
      allStories$call = c('germline', as.character(snpStories$call), as.character(cnvStories$call))
      allStories$stories = as.matrix(rbind(matrix(rep(1, length(qs)), nrow=1), snpStories$stories, cnvStories$stories))
      allStories$errors = as.matrix(rbind(matrix(rep(0, length(qs)), nrow=1), snpStories$errors, cnvStories$errors))
      rownames(allStories$stories) = rownames(allStories$errors) = rownames(allStories)
      colnames(allStories$stories) = colnames(allStories$errors) = ts

      #clusters stories into clones
      clusteredStories = storiesToCloneStories(allStories, minDistance=cloneDistanceCut, cpus=cpus, manualStoryMerge=manualStoryMerge, variants=variants)
      germlineCluster = which(apply(clusteredStories$cloneStories$stories + 1e-3 > 1, 1, all) &
        apply(clusteredStories$cloneStories$errors-1e-5 < 0, 1, all))
      rownames(clusteredStories$cloneStories)[germlineCluster] = clusteredStories$cloneStories$call[germlineCluster] = 'germline'
      rownames(clusteredStories$cloneStories$stories)[germlineCluster] = rownames(clusteredStories$cloneStories$errors)[germlineCluster] = 'germline'
      names(clusteredStories$storyList)[germlineCluster] = 'germline'
      catLog('got', nrow(clusteredStories$cloneStories), 'stories...')

      #combine the clustered stories into clonal evolution (figure out which is subclone of which)
      catLog('decide subclone structure..')
      cloneTree = findCloneTree(clusteredStories$cloneStories)

      catLog('removing inconsistent dodgy clones..')
      consistentTree = makeTreeConsistent(cloneTree, clusteredStories$cloneStories, clusteredStories$storyList)
      consistentClusteredStories = clusteredStories
      consistentClusteredStories$cloneTree = consistentTree
      consistentClusteredStories$cloneStories = consistentClusteredStories$cloneStories[clonesInTree(consistentTree),]
      consistentClusteredStories$storyList = consistentClusteredStories$storyList[clonesInTree(consistentTree)]
      catLog('got', nrow(consistentClusteredStories$cloneStories), 'consistent stories...')

      
      #add in previously filtered SNV and CNA stories if they fit with the found clones
      #also reassign anchor mutations to best fitting clone.
      #accept somatic SNVs down to 0.5 somaticP for these
      if ( any(!normals[ts]) ) {
        somaticMx = do.call(cbind, lapply(qs[!normals[ts]], function(q) q$somaticP > 0.5))
        somatic = apply(somaticMx, 1, any)
        catLog('Found ', sum(somatic), ' somatic SNVs to track. These will be linked to the found clones.\n', sep='')
      }
      else
        somatic = rep(FALSE, nrow(qs[[1]]))
      if ( any(normals[ts]) & rareGermline ) {
        rareGermlineMx = do.call(cbind, lapply(qs[normals[ts]], function(q) q$somaticP > 0.5 & q$severity < 10 & !is.na(q$severity)))
        isRareGermline = apply(rareGermlineMx, 1, any)
        catLog('Found ', sum(isRareGermline), ' rare germline variants to track.\n', sep='')
      }
      else
        isRareGermline = rep(FALSE, nrow(qs[[1]]))
      somaticQs = lapply(qs, function(q) q[isRareGermline | somatic,])
      rareGermlineNames = rownames(qs[[1]])[isRareGermline]
      
      snpStories = findSNPstories(somaticQs, cnvs[ts], normals[ts], filter=F, germlineVariants=rareGermlineNames)
      cnvStories = getCNVstories(cnvs[ts], normals[ts], genome, filter=F)
      allStories = data.frame(row.names=c('germline', rownames(snpStories), rownames(cnvStories)), stringsAsFactors=F)
      allStories$x1 = c(NA, snpStories$x1, cnvStories$x1)
      allStories$x2 = c(NA, snpStories$x2, cnvStories$x2)
      allStories$call = c('germline', as.character(snpStories$call), as.character(cnvStories$call))
      allStories$stories = as.matrix(rbind(matrix(rep(1, length(qs)), nrow=1), snpStories$stories, cnvStories$stories))
      allStories$errors = as.matrix(rbind(matrix(rep(0, length(qs)), nrow=1), snpStories$errors, cnvStories$errors))
      rownames(allStories$stories) = rownames(allStories$errors) = rownames(allStories)
      colnames(allStories$stories) = colnames(allStories$errors) = ts

      #reassign mutations to the consistent clones
      consistentClusteredStories$storyList = lapply(consistentClusteredStories$storyList, function(l) c())
      consistentClusteredStories = mergeStories(consistentClusteredStories, allStories, germlineVariants=rareGermlineNames)

      #redo clustering of mutations, this time using unfiltered mutations, and no consistency contstraints.
      use = rownames(allStories)
      if ( length(use) > maxStories ) {
        catLog('restricting to ', maxStories, ' random mutations in dodgy clustering the interest of time.\n', sep='')
        use = use[c(1,sample(2:length(use), maxStories))]
      }
      clusteredStories = storiesToCloneStories(allStories[use,], minDistance=cloneDistanceCut, cpus=cpus)
      germlineCluster = which(apply(clusteredStories$cloneStories$stories + 1e-3 > 1, 1, all) &
        apply(clusteredStories$cloneStories$errors-1e-5 < 0, 1, all))
      rownames(clusteredStories$cloneStories)[germlineCluster] = clusteredStories$cloneStories$call[germlineCluster] = 'germline'
      rownames(clusteredStories$cloneStories$stories)[germlineCluster] = rownames(clusteredStories$cloneStories$errors)[germlineCluster] = 'germline'
      names(clusteredStories$storyList)[germlineCluster] = 'germline'
      cloneTree = findCloneTree(clusteredStories$cloneStories)


      #separate stories that fitted to a consistent clone
      allConsistentStories = allStories[unlist(consistentClusteredStories$storyList),]
      
      #pick out the variants that behave like germline, ie present clonaly in all samples
      if ( nrow(clusteredStories$cloneStories) > 1 ) {
        germlineVariants = consistentClusteredStories$storyList[[which(rownames(consistentClusteredStories$cloneStories) == 'germline')]]
        germlineVariants = germlineVariants[germlineVariants != 'germline']
      }
      else germlineVariants = c()

      consistentClusteredStories = renameClones(consistentClusteredStories)
      
      stories[[name]] = list('allConsistent'=allConsistentStories, 'consistentClusters'=consistentClusteredStories, 'germlineVariants'=germlineVariants, 'all'=allStories, 'clusters'=clusteredStories, 'cloneTree'=cloneTree, anchorStories=list('anchorSNVs'=anchorSNVs, 'anchorCNAs'=anchorCNAs))
      catLog('done!\n')
    }
  }

  #flag variants that ended up in the germline clone
  #If rareGermline is switched off, instead set somaticP to 0.
  catLog('\nMarking SNV that behave like germline SNPs..')
  variants$variants = lapply(variants$variants, function(q) {
    q$germline = rep(NA, nrow(q))
    return(q)
  })
  if ( length(timeSeries) > 0 ) {
    for ( ind in names(stories) ) {
      story = stories[[ind]]
      if ( length(story$germlineVariants) > 0 ) {
        SNVs = story$germlineVariants[grep('^[0-9]', story$germlineVariants)]
        if ( length(SNVs) > 0 ) {
          samples = timeSeries[[ind]]
          variants$variants[samples] = lapply(variants$variants[samples], function(q) {
            if ( rareGermline ) q$germline = rownames(q) %in% SNVs
            else q$somaticP[rownames(q) %in% SNVs] = 0
            return(q)
          })
        }
      }
    }
  }
  catLog('saving variants back to file...')
  allVariantSaveFile = paste0(Rdirectory, '/allVariants.Rdata')
  load(allVariantSaveFile)
  allVariants$variants = variants
  save('allVariants', file=allVariantSaveFile)
  catLog('done!\n')

  #remove all variants in the germline clone from the clonal tracking.
  #This is isn't done earlier so that the somaticP can be switched off first in above chunk.
  if ( !rareGermline ) {
  	catLog('Removing potential germline variants from the output...')
  	stories = lapply(stories, function(story) {
  		#extract and remove all germline mutations from the consistent analysis
  		#keep only the first mutation that is the artificial "germline" mutation used as bait.
  		story$consistentClusters$storyList$germline = 'germline'
  		story$allConsistent = story$allConsistent[
  		           !(rownames(story$allConsistent) %in% story$germlineVariants),]
  		story$anchorStories$anchorSNVs = story$anchorStories$anchorSNVs[
  		           !(rownames(story$anchorStories$anchorSNVs) %in% story$germlineVariants),]
  		story$anchorStories$anchorCNAs = story$anchorStories$anchorCNAs[
  		           !(rownames(story$anchorStories$anchorCNAs) %in% story$germlineVariants),]
  		story$germlineVariants = 'germline'
  		#repeat for dodgy analysis
  		dodgyGLmuts = story$clusters$storyList$germline
  		dodgyGLmuts = dodgyGLmuts[dodgyGLmuts != 'germline']
  		story$clusters$storyList$germline = 'germline'
		story$all = story$all[!(rownames(story$all) %in% dodgyGLmuts),]
		return(story)
  	})
    catLog('done!\n')
  }

  #return clustered and raw stories
  stories = list('stories'=stories, 'variants'=variants, 'normalVariants'=allVariants$normalVariants)
  catLog('Saving stories..')
  save('stories', file=saveFile)
  catLog('done!\n\n\n')
  return(stories)
}

findSNPstories = function(somaticQs, cnvs, normal, filter=T, germlineVariants=c()) {
  if ( filter ) {
    cov10 = rowMeans(do.call(cbind, lapply(somaticQs, function(q) q$cov))) >= 10
    somaticQs = lapply(somaticQs, function(q) q[cov10,])
  }
  if ( nrow(somaticQs[[1]]) == 0 ) {
    emptyMx = matrix(,nrow=0, ncol=length(somaticQs))
    ret = data.frame(x1=numeric(), x2=numeric(), call=character(), stories=emptyMx, errors=emptyMx, stringsAsFactors=F)
    return(ret)
  }
  somaticQs = superFreq:::findSNPclonalities(somaticQs, cnvs)
  clonality = matrix(sapply(somaticQs, function(q) q$clonality), ncol=length(somaticQs))
  clonalityError = matrix(sapply(somaticQs, function(q) q$clonalityError), ncol=length(somaticQs))
  ret = data.frame(x1=somaticQs[[1]]$x, x2=somaticQs[[1]]$x, call=rownames(somaticQs[[1]]), row.names=rownames(somaticQs[[1]]), stringsAsFactors=F)
  clonality[which(is.na(clonality), arr.ind=T)] = 0
  clonalityError[which(is.na(clonalityError) | is.infinite(clonalityError), arr.ind=T)] = 10
  ret$stories = clonality
  ret$errors = clonalityError
  rownames(ret$stories) = rownames(ret$errors) = rownames(ret)
  
  colnames(ret$stories) = colnames(ret$errors) = names(somaticQs)
  if ( !filter ) {
    if ( any(normal) & nrow(ret) > 0 ) {
      allowed = rownames(ret) %in% germlineVariants
      presentInNormal = ret$stories[,normal,drop=F] > ret$errors[,normal,drop=F] | ret$stories[,normal,drop=F] > 0.2
      presentInNormal = apply(presentInNormal, 1, any)
      if ( class(presentInNormal) == 'matrix' ) presentInNormal = apply(presentInNormal, 1, any)
      ret = ret[!presentInNormal | allowed,,drop=F]
      catLog('Filtered ', sum(presentInNormal & !allowed), ' present in normal stories.\n', sep='')
    }
    return(ret)
  }
  
  allSmall = rowSums(is.na(ret$errors) | ret$stories - ret$errors*2 < 0 | ret$errors > 0.2) == ncol(ret$stories)
  ret = ret[!allSmall,,drop=F]
  uncertain = rowMeans(ret$errors) > 0.2
  ret = ret[!uncertain,,drop=F]
  indel = grepl('[-\\+]', rownames(ret))
  ret = ret[!indel,,drop=F]
  if ( any(normal) & nrow(ret) > 0 ) {
    presentInNormal = ret$stories[,normal,drop=F] > ret$errors[,normal,drop=F] | ret$stories[,normal,drop=F] > 0.2
    presentInNormal = apply(presentInNormal, 1, any)
    if ( class(presentInNormal) == 'matrix' ) presentInNormal = apply(presentInNormal, 1, any)
    ret = ret[!presentInNormal,,drop=F]
    catLog('Filtered ', sum(allSmall), ' small, ', sum(uncertain), ' uncertain, ', sum(indel), ' indel and ', sum(presentInNormal), ' present in normal stories.\n', sep='')
  }
  else catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain and ', sum(indel), ' indel stories.\n', sep='')

  return(ret)
}

findLocalCNV = function(qs, cnvs) {
  x = qs[[1]]$x
  for ( i in 1:length(cnvs) ) {
    call = sapply(x, function(X) c(cnvs[[i]]$clusters$call[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 'AB')[1])
    clonality = sapply(x, function(X) c(cnvs[[i]]$clusters$clonality[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 1)[1])
    clonalityError = sapply(x, function(X) c(cnvs[[i]]$clusters$clonalityError[which(X > cnvs[[i]]$clusters$x1-300 & X < cnvs[[i]]$clusters$x2+300)], 0)[1])
    if ( any(grepl('CL', call)) ) {
      call[grepl('CL', call)] = 'AB'
      clonality[grepl('CL', call)] = 1 - clonality[grepl('CL', call)]
    }
    qs[[i]]$CNV = call
    qs[[i]]$CNVclonality = clonality
    qs[[i]]$CNVclonalityError = clonalityError
  }
  return(qs)
}

#takes a quality variant object, that has gone through findLocalCNV, and adds a clonality and clonalityError column
findSNPclonalities = function(somaticQs, cnvs) {
  if ( nrow(somaticQs[[1]]) == 0 ) return(somaticQs)
  somaticQs = findLocalCNV(somaticQs, cnvs)
  somaticQs = lapply(somaticQs, function(q) {
    nA = nchar(q$CNV) - nchar(gsub('A', '', q$CNV))
    nB = nchar(q$CNV) - nchar(gsub('B', '', q$CNV))

    noCNV = q$CNV %in% c('AB', 'AB?', 'AB??')

    #calculate expected frequencies range if the SNP is on the A or B allele of the CNV, or in the AB background.
    cloneMin = ifelse(noCNV, 0, pmax(0, q$CNVclonality-q$CNVclonalityError*1.5))
    cloneMax = ifelse(noCNV, 0, pmin(1, q$CNVclonality+q$CNVclonalityError*1.5))
    fAmin = nA*cloneMin/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fAmax = nA*cloneMax/((nA+nB)*cloneMax + 2*(1-cloneMax))
    fBmin = nB*cloneMin/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fBmax = nB*cloneMax/((nA+nB)*cloneMax + 2*(1-cloneMax))
    fNmax = (1-cloneMin)/((nA+nB)*cloneMin + 2*(1-cloneMin))
    fNmin = (1-cloneMax)/((nA+nB)*cloneMax + 2*(1-cloneMax))
    
    #calculate p-values for the SNV living on the different alleles.
    f = superFreq:::refUnbias(q$var/q$cov)
    f[q$cov==0] = 0
    closestFA = ifelse(f < fAmin, fAmin, ifelse(f > fAmax, fAmax, f))
    closestFB = ifelse(f < fBmin, fBmin, ifelse(f > fBmax, fBmax, f))
    closestFN = ifelse(f < fNmax, f, fNmax)
    pA = superFreq:::pBinom(q$cov, q$var, closestFA)
    pB = superFreq:::pBinom(q$cov, q$var, closestFB)
    pN = superFreq:::pBinom(q$cov, q$var, closestFN)

    clonalityA = pmin(1, ifelse(nA==0 & f==0, 0, 2*f/(nA+f*(2-nA-nB))))
    clonalityB = pmin(1, ifelse(nB==0 & f==0, 0, 2*f/(nB+f*(2-nA-nB))))
    clonalityN = pmin(1, 1 - (1 - f*2)/(f*(nA+nB-2)+1))

    #frequency error estimate: add up the poissonian width with the RIB as independent normal error sources.
    fErr = sqrt(superFreq:::frequencyError(q$var, q$cov, p0=0.10)^2 + q$RIB^2)
    fErr[q$cov==0] = 1
    #propagate to clonality
    clonalityHighA = 2/(2+nA/(f+fErr)-nA-nB)
    clonalityHighB = 2/(2+nB/(f+fErr)-nA-nB)
    clonalityHighN = 1 - (1 - (f+fErr)*2)/(pmin(1, f+fErr)*(nA+nB-2)+1)
    #Not allowing the lower limit to go into negatives reduces the error estimate of low frequencies
    #and is associated with a prior bias towards clonality 0 for low frequencies.
    clonalityLowA = pmax(0, 2/(2+nA/(f-fErr)-nA-nB))
    clonalityLowB = pmax(0, 2/(2+nB/(f-fErr)-nA-nB))
    clonalityLowN = pmax(0, 1 - (1 - (f-fErr)*2)/((f-fErr)*(nA+nB-2)+1))
      
    clonalityErrorA = abs(clonalityHighA - clonalityLowA)/2
    clonalityErrorB = abs(clonalityHighB - clonalityLowB)/2
    clonalityErrorN = abs(clonalityHighN - clonalityLowN)/2
    #clonalityError[f==0] = fErr[f==0]

    #if the B allele is lost, and the SNV is clonal in the N-cells, assume that the SNV
    #was present in the cells with the lost B-allele as well, giving a clonality of 1.
    lostSNV = nB==0 & abs(clonalityN+q$CNVclonality - 1) < q$CNVclonalityError
    if ( any(lostSNV) ) clonalityErrorN[lostSNV] = 1

    consistentA = pA > 0.05 & clonalityLowA <= cloneMax
    consistentB = pB > 0.05 & clonalityLowB <= cloneMax
    consistentN = pN > 0.05 & clonalityLowN <= 1-cloneMin
    
    #if no allele fits with the point mutation and CNA present in mostly the same cells
    #then check for SNV present in subclone of CNAs with a single allele
    #inconsistent = !consistentA & !consistentB & !consistentN
    #if ( any(inconsistent) ) {
	#	consistentSubA = f < fAmax & nA == 1
	#	consistentSubB = f < fBmax & nB == 1
	#	consistentSubN = f < fNmax
	#	consistentA[inconsistent] = consistentSubA[inconsistent]
	#	consistentB[inconsistent] = consistentSubB[inconsistent]
	#	consistentN[inconsistent] = consistentSubN[inconsistent]
    #}
    
    homeAllele =
      ifelse(q$CNV %in% c('AB', 'AB?', 'AB??', 'CL'),
             'N',
             ifelse(consistentA & (!consistentB | clonalityB <= clonalityA) & (!consistentN | clonalityN <= clonalityA),
                    'A',
                    ifelse(consistentB & (!consistentN | clonalityN <= clonalityB),
                           'B',
                           'N')))
    clonality = ifelse(homeAllele == 'A', clonalityA, ifelse(homeAllele == 'B', clonalityB, clonalityN))
    homeAllele[f==0] = 'N'
    clonality[f==0] = 0
    clonalityError = ifelse(homeAllele == 'A', clonalityErrorA, ifelse(homeAllele == 'B', clonalityErrorB, clonalityErrorN))

    #check for calls where several different home alleles are possible, and increase
    #error estimate accordingly
    uncertainCalls = which(consistentA + consistentB + consistentN > 1)
    if ( length(uncertainCalls) > 0 ) {
      clonalityVar = sapply(uncertainCalls, function(i) var(c(clonalityA[i], clonalityB[i], clonalityN[i])[c(consistentA[i], consistentB[i], consistentN[i])]))
      clonalityError[uncertainCalls] = sqrt(clonalityError[uncertainCalls]^2 + clonalityVar)
    }

    #check if the snv can be a germline SNP, ie present both on the N and one of the A or B alleles.
    #if so, make sure that the clonality error reaches 100%.
    SNPfA = (q$CNVclonality*nA + 1-q$CNVclonality)/(2*(1-q$CNVclonality) + q$CNVclonality*(nA+nB))
    SNPfB = (q$CNVclonality*nB + 1-q$CNVclonality)/(2*(1-q$CNVclonality) + q$CNVclonality*(nA+nB))
    canBeSNPA = superFreq:::pBinom(q$cov, q$var, SNPfA)
    canBeSNPB = superFreq:::pBinom(q$cov, q$var, SNPfB)
    canBeSNP = canBeSNPA > 0.1 | canBeSNPB > 0.1
    clonalityError[canBeSNP] = pmax(clonalityError[canBeSNP], (1-clonality)[canBeSNP])


    q$homeAllele = homeAllele
    q$clonality = superFreq:::noneg(clonality)
    q$clonalityError = abs(clonalityError)
    return(q)
    })
  return(somaticQs)
}

frequencyError = function(var, cov, p0=0.15, reportBothEnds=F) {
  covOri = cov
  varOri = var

  converged = cov == 0
  bsolutions = ifelse(cov==0, 1, 0)
  var = var[!converged]
  cov = cov[!converged]
  b = var - round(sqrt(var))
  p = ifelse(b < 0, 0, 1-pbinom(var-1, cov, noneg(b)/cov))
  while ( any(!converged) ) {
    tooLow = p < p0
    bnew = b + sign(p0-p)
    pnew = ifelse(bnew < 0, 0, 1-pbinom(var-1, cov, noneg(bnew)/cov))
    solved = sign(p-p0) != sign(pnew-p0) 
    #handle converged cases
    if ( any(solved) ) {
      bsolutions[which(!converged)[solved]] = ((b*(pnew-p0)+bnew*(p0-p))/(pnew-p))[solved]
      converged[which(!converged)[solved]] = rep(T, sum(solved))
    }
    #update not converged cases, removing solved values
    p = pnew[!solved]
    b = bnew[!solved]
    cov = cov[!solved]
    var = var[!solved]
  }

  var = varOri
  cov = covOri
  converged = cov == 0
  tsolutions = ifelse(cov==0, 1, 0)
  var = var[!converged]
  cov = cov[!converged]
  t = var + round(sqrt(var))
  p = ifelse(t > cov, 0, pbinom(var, cov, pmin(cov, t)/cov))
  while ( any(!converged) ) {
    tooLow = p < p0
    tnew = t + sign(p-p0)
    pnew = ifelse(tnew > cov, 0, pbinom(var, cov, pmin(cov, tnew)/cov))
    solved = sign(p-p0) != sign(pnew-p0) 
    #handle converged cases
    if ( any(solved) ) {
      tsolutions[which(!converged)[solved]] = ((t*(pnew-p0)+tnew*(p0-p))/(pnew-p))[solved]
      converged[which(!converged)[solved]] = rep(T, sum(solved))
    }
    #update not converged cases, removing solved values
    p = pnew[!solved]
    t = tnew[!solved]
    cov = cov[!solved]
    var = var[!solved]
  }

  if ( reportBothEnds ) return(cbind(pmax(0,bsolutions/covOri), pmin(1,tsolutions/covOri)))
  error = pmax(tsolutions-varOri, varOri-bsolutions)/covOri
  return(error)
}

getCNVstories = function(cnvs, normal, genome, filter=T) {
  regions = splitRegions(cnvs)
  catLog(nrow(regions), ' regions..', sep='')
  events = splitEvents(cnvs, regions)
  catLog(nrow(regions), ' events.\nExtracting stories, comparing SNP shift directions..', sep='')
  stories = cnvsToStories(cnvs, events, normal, genome, filter=filter)
  return(stories)
}

splitRegions = function(cnvs) {
  regions = data.frame('x1' = c(), 'x2' = c(), stringsAsFactors=F)
  cnvX1 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x1)))
  cnvX2 = unique(unlist(sapply(cnvs, function(cnv) cnv$cluster$x2)))
  x1 = x2 = -Inf
  while (x1 < max(cnvX1)) {
    x1 = min(cnvX1[cnvX1 >= x2])
    x2 = min(cnvX2[cnvX2 > x1])
    regions = rbind(regions, data.frame('x1'=x1, 'x2'=x2, stringsAsFactors=F))
  }
  return(regions)
}

splitEvents = function(cnvs, regions) {
  events = data.frame('x1' = c(), 'x2' = c(), 'call' = c(), stringsAsFactors=F)
  cnvX1 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x1))
  cnvX2 = unlist(sapply(cnvs, function(cnv) cnv$cluster$x2))
  cnvCalls = unlist(sapply(cnvs, function(cnv) cnv$cluster$call))
  for ( i in 1:nrow(regions) ) {
    calls = unique(cnvCalls[cnvX1 <= regions$x1[i] & cnvX2 >= regions$x2[i]])
    calls = calls[!(calls %in% c('AB', 'AB?', 'AB??'))]
    for ( call in calls ) events = rbind(events, data.frame('x1' = regions$x1[i], 'x2' = regions$x2[i], 'call' = call, stringsAsFactors=F))
  }
  events = events[!grepl('\\?', events$call),]
  return(events)
}

cnvsToStories = function(cnvs, events, normal, genome, filter=T) {
  stories = errors = sigmas = data.frame(stringsAsFactors=F)
  if ( nrow(events) > 0 ) {
    i=1
    while ( i <= nrow(events) ) {
      call = as.character(events$call[i])
      clonalities = extractClonalities(cnvs, events[i,])
      stories = rbind(stories, noneg(clonalities$clonality))
      errors = rbind(errors, clonalities$clonalityError)
      sigmas = rbind(sigmas, clonalities$sigma)
      if ( any(clonalities$clonality < 0) & any(clonalities$clonality > 0) & !any(is.na(clonalities$clonality)) ) {
        negativeClonalities = clonalities
        negativeClonalities$clonality = noneg(-negativeClonalities$clonality)
        negEvent = events[i,]
        negEvent$call = reverseCall(negEvent$call)
        events = rbind(events, negEvent)[order(c(1:nrow(events), i+0.5)),]
        stories = rbind(stories, noneg(negativeClonalities$clonality))
        errors = rbind(errors, negativeClonalities$clonalityError)
        sigmas = rbind(sigmas, negativeClonalities$sigma)
        i = i + 1
      }
      i = i + 1
    }

    if ( nrow(events) > 1 ) {
      Nevents = nrow(events)
      keep = rep(T, Nevents)
      for ( x1 in unique(events$x1) ) {
        rows = events$x1 == x1
        if ( sum(rows) == 0 ) next
        subStories = stories[rows,,drop=F]
        subSigmas = sigmas[rows,,drop=F]
        subKeep = keep[rows]
        while ( !any(is.na(colSums(subStories[subKeep,,drop=F]))) & any(colSums(subStories[subKeep,,drop=F]) > 1) ) {
          presence = rowMeans(subStories[subKeep,,drop=F])
          inconsistency = sqrt(rowMeans(subSigmas[subKeep,,drop=F]^2))
          score = presence - inconsistency/3
          subKeep[subKeep][which.min(score)] = F
        }
        keep[rows] = subKeep
      }
      catLog(sum(!keep), 'overlapping stories filtered..')
      events = events[keep,,drop=F]
      stories = stories[keep,,drop=F]
      errors = errors[keep,,drop=F]
      sigmas = sigmas[keep,,drop=F]
    }
  }
  else {
    events$stories=matrix(,nrow=0, ncol=length(cnvs))
    events$errors=matrix(,nrow=0,ncol=length(cnvs))
    return(events)
  }
  colnames(errors) = colnames(stories) = names(cnvs)
  ret = events
  rownames(ret) = make.names(paste0('chr', xToChr(ret$x1, genome), '-', events$call), unique=T)
  ret$stories = as.matrix(stories)
  ret$errors = as.matrix(errors)
  ret$errors[is.na(ret$errors)] = 1

  if ( !filter ) return(ret)
  
  allSmall = rowSums(ret$stories - ret$errors*1.5 < 0.15 | ret$errors > 0.2) == ncol(ret$stories)
  ret = ret[!allSmall,,drop=F]
  uncertain = rowMeans(ret$errors) > 0.2
  ret = ret[!uncertain,,drop=F]
  notSignificant = rowSums(abs(ret$stories) < ret$errors*1.5 | is.na(ret$errors)) >= ncol(ret$stories)
  ret = ret[!notSignificant,,drop=F]
  smallRegion = ret$x2 - ret$x1 < 5e6 & !(ret$call == 'CL' & apply(ret$stories, 1, max) > 0.7)
  ret = ret[!smallRegion,,drop=F]
  if ( any(normal) & nrow(ret) > 0 ) {
    presentInNormal = ret$stories[,normal,drop=F] > ret$errors[,normal,drop=F] | ret$stories[,normal,drop=F] > 0.2
    presentInNormal = apply(presentInNormal, 1, any)
    ret = ret[!presentInNormal,,drop=F]
    catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain, ', sum(smallRegion), ' small region and ', sum(presentInNormal), ' present in normal stories.\n', sep='')
  }
  else catLog('Filtered ', sum(allSmall) , ' small, ', sum(uncertain), ' uncertain and ', sum(smallRegion), ' small region stories.\n', sep='')
  falseSNPcalls = ret$call == 'AA' & (ret$x2 - ret$x1 < 2e6 | rowMeans(ret$errors) > 0.15 )
  ret = ret[!falseSNPcalls,,drop=F]
  catLog('Filtered ', sum(falseSNPcalls) , ' stories that are potentially based on false SNPs.\n', sep='')
  return(ret)
}

extractClonalities = function(cnvs, event) {
  nSample = length(cnvs)
  subFreqs = lapply(cnvs, function(cs) cs$eFreqs[cs$eFreqs$x >= event$x1 & cs$eFreqs$x <= event$x2,])
  subFreqsMirror = lapply(cnvs, function(cs) {
    ret=cs$eFreqs[cs$eFreqs$x >= event$x1 & cs$eFreqs$x <= event$x2,]
    ret$var = mirrorDown(ret$var, ret$cov)
    return(ret)
    })
  subCR = do.call(rbind, lapply(1:length(cnvs), function(i) {
    cs = cnvs[[i]]
    efs = subFreqsMirror[[i]]
    ret = mergeToOneRegion(cs$CR[cs$CR$x1 >= event$x1 & cs$CR$x2 <= event$x2,,drop=F], efs)
  }))
  fM = callTofM(as.character(event$call))
  prior = callPrior(as.character(event$call))
  ret = as.data.frame(do.call(rbind, lapply(1:nSample, function(sample) {
    efs=subFreqsMirror[[sample]]
    iscnv = isCNV(cluster=subCR[sample,], efs=efs, M=fM[2], f=fM[1], prior=prior)
    #theoreticalError = getTheoreticalError(cluster=subCR[sample,], call=as.character(event$call), efs=efs)
    #iscnv$clonalityError = sqrt(iscnv$clonalityError^2 + theoreticalError^2)
    bestSigma = getBestSigma(cluster=subCR[sample,], call=as.character(event$call), efs=efs)
    return(c(iscnv$call, iscnv$clonality, iscnv$sigma, iscnv$clonalityError, bestSigma))
  })), stringsAsFactors=F)
  names(ret) = c('call', 'clonality', 'sigma', 'clonalityError', 'bestSigma')
  direction = rep(0, nSample)

  if ( fM[1] != 0.5 ) {
    freqX = unique(unlist(lapply(subFreqs, function(freq) freq$x)))
    if ( length(freqX) > 0 ) {
      directionScores = do.call(cbind, lapply(1:nSample, function(i) freqToDirectionProb(subFreqs[[i]], freqX, fM, ret$clonality[i])))
      strongestSample = which(colsums(directionScores^2) >= 0.999999*max(colsums(directionScores^2)))[1]
      direction = sapply(1:nSample, function(i) {
        sN = scalarNorm(directionScores[,strongestSample], directionScores[,i])
        power = sqrt(sum(directionScores[,i]^2))
        if ( (is.na(power) == TRUE) || power < 10 ) return(NA)  #if not enough power to tell up or down, nevermind.
        return(sN)
      })
    }
  }

  #if there are other much better calls, this call probably isnt present
  wrongCall = ret$sigma > 5 & ret$bestSigma < 3
  ret$clonality[wrongCall] = 0

  #if no call gives a good fit, just up the uncertainty depending on how bad fit.
  #messedUpCall = ret$sigma > 5 & !wrongCall
  #ret$clonalityError[messedUpCall] = sqrt(ret$clonalityError[messedUpCall]^2 + ((ret$sigma[messedUpCall]-5)/5)^2)
  ret$clonalityError = sqrt(ret$clonalityError^2 + ((ret$sigma - 2)/20)^2)
  
  ret$clonality[!is.na(direction) & direction < 0] = -ret$clonality[!is.na(direction) & direction < 0]
  if ( !any(ret$clonality > 0) & any(ret$clonality < 0) ) ret$clonality = - ret$clonality
  
  return(ret)
}

mergeToOneRegion = function(cR, eFreqs) {
  if ( nrow(cR) == 0 )
    return(data.frame(x1=NA, x2=NA, M=0, width=1000, df=10, var=0, cov=0,
                      pHet=0.5, pAlt=0.5, odsHet=0.5, f=NA, ferr=NA))
  cR$x1[1] = min(cR$x1)
  cR$x2[1] = max(cR$x2)
  efs = eFreqs[eFreqs$x > cR$x1[1] & eFreqs$x < cR$x2[1],]
  cR$var[1] = sum(cR$var)
  cR$cov[1] = sum(cR$cov)
  cR$M[1] = sum(cR$M/cR$width^2)/sum(1/cR$width^2)
  cR$width[1] = 1/sqrt(sum(1/cR$width^2))
  cR$Nsnps[1] = sum(cR$Nsnps)
  cR = cR[1,]
  cf = correctedFrequency(cR, efs)
  cR$f = cf[,'f']
  cR$ferr = cf[,'ferr']
  return(cR)
}

freqToDirectionProb = function(freq, freqX, fM, clonality) {
  f = fM[1]
  M = fM[2]
  if ( f == 0.5 ) return(rep(0, length(freqX)))
  
  upProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    fClone = (clonality*(2*f*2^M-1) + 1)/(clonality*(2*2^M - 2) + 2)
    return(pBinom(freq$cov[j], freq$var[j], 1-fClone))
  }
  downProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    fClone = (clonality*(2*f*2^M-1) + 1)/(clonality*(2*2^M - 2) + 2)
    return(pBinom(freq$cov[j], freq$var[j], fClone))
  }
  nullProb = function(x) {
    j = which(freq$x == x)
    if ( length(j) == 0 ) return(1)
    return(pBinom(freq$cov[j], freq$var[j], 0.5))
  }
  u = 1e-5 + sapply(freqX, function(x) upProb(x))
  d = 1e-5 + sapply(freqX, function(x) downProb(x))
  n = 1e-5 + sapply(freqX, function(x) nullProb(x))
  return((1-n)*log(u/d))
}



#takes a dataframe of stories and groups them into subclone stories. returns a data frame of the subclone stories
#and a list of dataframes for the individual stories in each subclone.
storiesToCloneStories = function(stories, storyList='',
  minDistance=-qnorm(0.01), cpus=1, manualStoryMerge=F, variants=NA) {
  if ( (storyList ==  '')[1] ) storyList = as.list(rownames(stories))
  if ( length(storyList) < 2 )
    return(list('cloneStories'=stories, 'storyList'=storyList))

  #if a lot of mutations, merge them in batches, as the algorithm scales as O(N^2) otherwise
  #group less in this first pass (minDistance*0.5), and then group as specified last round.
  batchSize = round(pmin(1000, pmax(100, 1e6/length(storyList))))
  if ( length(storyList) > batchSize ) catLog('Cluster mutations by batch. Remaining clusters: ', sep='')
  while ( length(storyList) > batchSize ) {
    catLog(length(storyList), '..', sep='')
    first300 = storiesToCloneStories(stories=stories, storyList=storyList[1:batchSize],
      minDistance=minDistance*0.5, cpu=cpus)
    storyList = c(first300$storyList, storyList[(batchSize+1):length(storyList)])
    batchSize = round(pmin(1000, pmax(100, 1e6/length(storyList))))
  }

  summariseClusters = function() {
    st = do.call(rbind, lapply(storyList, function(rows) {
      err = stories$errors[rows,,drop=F]
      if ( any(err<=0) ) err[err<=0] = rep(min(c(1,err[err>0]))/1e6, sum(err<=0))
      st = stories$stories[rows,,drop=F]
      w = t(1/t(err^2)/colsums(1/err^2))
      mean = colsums(st*w)
      ret = matrix(mean, nrow=1)
    }))
    err = do.call(rbind, lapply(storyList, function(rows) {
      err = stories$errors[rows,,drop=F]
      err = 1/sqrt(colsums(1/err^2))
      ret = matrix(err, nrow=1)
    }))
    if ( nrow(st) > 0 ) rownames(err) = rownames(st) = 1:length(storyList)
    colnames(err) = colnames(st) = colnames(stories$stories)
    clusters = data.frame('call'=rep('clone', length(storyList)), 'x1'=rep(NA, length(storyList)), 'x2'=rep(NA, length(storyList)), row.names=1:length(storyList), stringsAsFactors=F)
    clusters$stories = st
    clusters$errors = err
    return(clusters)
  }
  
  distance = do.call(rbind, mclapply(1:length(storyList), function(i) c(sapply(1:i, function(j) pairScore(stories, storyList[[i]], storyList[[j]])), rep(0, length(storyList)-i)), mc.cores=cpus))
  distance = distance + t(distance)
  distance = distance + ifelse(row(distance) == col(distance), (max(distance, minDistance)+1), 0)

  loops = 0
  while ( any(distance < minDistance) ) {
    merge = which(distance == min(distance), arr.ind=TRUE)[1,]

    if ( manualStoryMerge & loops == 0 ) {
      layout(matrix(1:2, nrow=2))
      cat('\n\nUSER INPUT:\nTop left: Merging story cluster ', merge[1], ' (blue) and ', merge[2], ' (red) at a distance of ', min(distance), '.\n', sep='')
      cat('There are ', nrow(distance), ' clusters.\n', sep='')
      i=merge[1];j=merge[2];plotStories(stories[c(storyList[[i]], storyList[[j]]),], variants, col=mcri(c(rep('blue', length(storyList[[i]])), rep('red', length(storyList[[j]])))), lty=1, main=paste0('next up to merge: ', merge[1], ' (blue) and ', merge[2], ' (red)'), setPar=F)
      clusters = summariseClusters()
      plotStories(clusters, variants, main = 'clusters right now', setPar=F)
      if ( nrow(distance) < 15 ) {
        cat('Clone distance matrix:\n')
        print(distance)
      }
      layout(1)

      a = as.numeric(readline('\nType number of merges to perform, or blank for this merge only.\nType 0 to not merge these two clusters.\nType -1 to reset and start over with clustering.\nType any number smaller than -1 to stop clustering and continue downstream analysis.\nInput:'))
      if ( is.na(a) ) a = 1
      if ( a < -1 ) {
        distance = Inf
        next
      }
      else if ( a == 0 ) {
        distance[merge[1], merge[2]] = distance[merge[2], merge[1]] = Inf
        next
      }
      else if ( a == -1 ) {
        storyList=as.list(rownames(stories))
        distance = do.call(rbind, mclapply(1:length(storyList), function(i) c(sapply(1:i, function(j) pairScore(stories, storyList[[i]], storyList[[j]])), rep(0, length(storyList)-i)), mc.cores=cpus))
        distance = distance + t(distance)
        distance = distance + ifelse(row(distance) == col(distance), (max(distance)+1), 0)
        next
      }
      else if ( a == round(a) & a > 0 ) loops = a - 1
      else loops = 0
    }
    else if ( manualStoryMerge ) loops = loops - 1
    
    storyList[[merge[1]]] = c(storyList[[merge[1]]], storyList[[merge[2]]])
    distance[merge[1],] = distance[,merge[1]] = sapply(1:length(storyList), function(j)
                                      pairScore(stories, storyList[[merge[1]]], storyList[[j]]) + ifelse(j==merge[1], (max(distance)+1), 0))
    distance = distance[-merge[2], -merge[2], drop=F]
    storyList = storyList[-merge[2]]
  }

  clusters = summariseClusters()
  names(storyList) = rownames(clusters)
  
  return(list('cloneStories'=clusters, 'storyList'=storyList))
}

#The metric on stories, used for clustering similar stories into subclones.
pairScore = function(stories, is, js) {
  if ( identical(is, js) ) return(0)
  #if ( length(is) == 1 ) rms1 = noneg(0.5 - mean(stories$errors[is,]))
  #else {
  #  err1 = stories$errors[is,]
  #  if ( any(err1<=0) ) err1[err1<=0] = rep(min(c(1,err1[err1>0]))/1e6, sum(err1<=0))
  #  st1 = stories$stories[is,]
  #  w1 = t(1/t(err1^2)/colsums(1/err1^2))
  #  mean1 = colSums(st1*w1)
  #  sigma1 = abs(t(t(st1)-mean1))/err1
  #  rms1 = max(sqrt(colmeans(sigma1^2)))
  #}
  #if ( length(js) == 1 ) rms2 = noneg(0.5 - mean(stories$errors[js,]))
  #else {
  #  err2 = stories$errors[js,]
  #  if ( any(err2<=0) ) err2[err2<=0] = rep(min(c(1,err2[err2>0]))/1e6, sum(err2<=0))
  #  st2 = stories$stories[js,]
  #  w2 = t(1/t(err2^2)/colsums(1/err2^2))
  #  mean2 = colSums(st2*w2)
  #  sigma2 = abs(t(t(st2)-mean2))/err2
  #  rms2 = max(sqrt(mean(sigma2^2)))
  #}
  #unpairedRms = sqrt(rms1^2+rms2^2)

  err = stories$errors[c(is,js),,drop=F]
  if ( any(err<=0) ) err[err<=0] = rep(min(c(1,err[err>0]))/1e6, sum(err<=0))
  if ( any(is.infinite(err)) ) err[is.infinite(err)] = rep(1e6, sum(is.infinite(err)))
  st = stories$stories[c(is,js),,drop=F]
  w = t(1/t(err^2)/colsums(1/err^2))
  mean = colSums(st*w)
  sigma = abs(t(t(st)-mean))/err
  #rms = max(sqrt(colmeans(sigma^2)))

  #p = min(apply(pnorm(-sigma)*2, 2, function(ps) min(p.adjust(ps, method='fdr'))))
  #p1 = min(apply(pnorm(-sigma[1:length(is),])*2, 2, function(ps) min(p.adjust(ps, method='fdr'))))
  #p2 = min(apply(pnorm(-sigma[(length(is)+1):nrow(sigma),])*2, 2, function(ps) min(p.adjust(ps, method='fdr'))))
  #pF = min(apply(pnorm(-sigma)*2, 2, function(ps) fisherTest(ps)['pVal']))
  pF1 = min(apply(pnorm(-sigma)*2, 2, function(ps) fisherTest(ps[1:length(is)])['pVal']))
  pF2 = min(apply(pnorm(-sigma)*2, 2, function(ps) fisherTest(ps[(length(is)+1):length(ps)])['pVal']))
  return(-qnorm(min(pF1, pF2)/2))
  
  #return(rms + noneg(rms - unpairedRms) )
}


#helper function that decided which subclones are subclones of each other.
findCloneTree = function(cloneStories, storyList) {
  #add the purity as a story (this may or may not already be present)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(cloneStories$stories[,i]))
  #check which clones are subclones of which others, allowing an error bar from each clone.
  subcloneMx = as.matrix(apply(abs(cloneStories$stories) + cloneStories$errors*1.5, 1,
    function(story1) apply(abs(cloneStories$stories) - cloneStories$errors*1.5, 1,
                           function(story2) all(story1 > story2))))
  greaterSum = as.matrix(apply(abs(cloneStories$stories), 1,
    function(story1) apply(abs(cloneStories$stories), 1,
                           function(story2) sum(story1) > sum(story2))))
  subcloneMx = subcloneMx & greaterSum
  subcloneMx = enforceTransitive(subcloneMx)
  subcloneStories = cloneStories$stories
  rownames(subcloneMx) = colnames(subcloneMx) = rownames(subcloneStories) = rownames(cloneStories)
  
  cloneTree = findChildren(subcloneMx, subcloneStories)

  return(cloneTree)
}
#helper function
findChildren = function(subcloneMx, subcloneStories) {
  #if no more subclones, return empty list
  if ( nrow(subcloneStories) == 0 ) return(list())

  #there will be at least one clone that is not a subclone (aka paranetless, the one with the largest clonality sum)
  cloneNames = rownames(subcloneStories)[which(rowSums(subcloneMx) == 0)]
  
  #score the parentless subclones from the sum of clonalities.
  cloneScores = rowSums(abs(subcloneStories))[cloneNames]
  cloneScores = sort(cloneScores, decreasing=T)

  #go through the parentless clones in order of score, each one recurring with its subclones.
  ret = list()
  for ( clone in names(cloneScores) ) {
    subclones = subcloneMx[,clone]
    if ( !any(subclones) ) ret[[clone]] = list()
    else ret[[clone]] = findChildren(subcloneMx[subclones, subclones, drop=F], subcloneStories[subclones,,drop=F])
    subcloneMx = subcloneMx[!subclones, !subclones, drop=F]
    subcloneStories = subcloneStories[!subclones,,drop=F]
  }

  return(ret)
}

#helper function, switching A and B to show that the other allele is affected.
reverseCall = function(call) {
  call = gsub('A', 'a', call)
  call = gsub('B', 'A', call)
  call = gsub('a', 'B', call)
  return(call)
}

#enforces transitivity on a subclone matrix by taking the transitive closure
enforceTransitive = function(mx) {
  for ( col in 1:ncol(mx) ) {
    mx[,col] = mx[,col] | rowsums(mx[,mx[,col],drop=F]) > 0
  }
  return(mx)
}

#add the filtered stories to the clone stories if consistent.
mergeStories = function(clusteredStories, filteredStories, germlineVariants=c()) {
  germlineOnly = rownames(filteredStories) %in% germlineVariants
  distance = sapply(1:nrow(clusteredStories$cloneStories), function(cloneRow) {
    isGermlineClone = rownames(clusteredStories$cloneStories)[cloneRow] == 'germline'
    sapply(1:nrow(filteredStories), function(storyRow) {
      ret = cloneStoryRMS(clusteredStories$cloneStories[cloneRow,], filteredStories[storyRow,])
      if ( germlineOnly[storyRow] & !isGermlineClone ) ret = Inf
      #be extra generous for rare germline variants matching the germline clone.
      if ( germlineOnly[storyRow] & isGermlineClone ) ret = ret/2
      return(ret)
    })
  })
  #if either as only one row, make it a 1-row or 1-column matrix.
  if (any(class(distance) != 'matrix'))
    distance = matrix(distance, nrow = nrow(filteredStories))
  closestDistance = apply(distance, 1, min)
  toMerge = which(closestDistance <= 1)
  bestClone = sapply(toMerge, function(row) which(distance[row,] == closestDistance[row])[1])
  for ( i in 1:length(toMerge) )
    clusteredStories$storyList[[bestClone[i]]] = c(clusteredStories$storyList[[bestClone[i]]], rownames(filteredStories)[toMerge[i]])

  return(clusteredStories)
}

#distance between a story and a clone
cloneStoryRMS = function(clone, story, systematicVariance=0.02) {
  errors = sqrt(clone$errors^2+story$errors^2)+0.02
  sigmas = (clone$stories - story$stories)/errors
  rms = sqrt(mean(sigmas^2))
  return(rms)
}

combineStories = function(stories1, stories2) {
  ret = data.frame(row.names=c(rownames(stories1), rownames(stories2)), stringsAsFactors=F)
  ret$x1 = c(stories1$x1, stories2$x1)
  ret$x2 = c(stories1$x2, stories2$x2)
  ret$call = c(as.character(stories1$call), as.character(stories2$call))
  ret$stories = as.matrix(rbind(stories1$stories, stories2$stories))
  ret$errors = as.matrix(rbind(stories1$errors, stories2$errors))
  rownames(ret$stories) = rownames(ret$errors) = rownames(ret)
  return(ret)
}


addTheoreticalCNVerrors = function(cnvs) {
  for ( sample in names(cnvs) ) {
    clusters = cnvs[[sample]]$clusters
    eFreqs = cnvs[[sample]]$eFreqs
    eFreqs$var = mirrorDown(eFreqs$var, eFreqs$cov)
    for ( row in 1:nrow(clusters) ) {
      call = clusters$call[row]
      if ( grepl('\\?', call) ) next
      if ( call == 'AB' ) next
      efs = eFreqs[eFreqs$x > clusters$x1[row] & eFreqs$x < clusters$x2[row],]

      theoreticalError = getTheoreticalError(clusters[row,], clusters$call[row], efs)
      clusters[row,]$clonalityError = sqrt(clusters[row,]$clonalityError^2 + theoreticalError^2)
    }

    cnvs[[sample]]$clusters = clusters
  }

  return(cnvs)
}

getTheoreticalError = function(cluster, call, efs) {
  isCall = isCNV(cluster, efs, callTofM(call)['M'],
    callTofM(call)['f'], callPrior(call), 5)
  bestSigma = isCall$sigma
  bestClonality = isCall$clonality
  
  secondBestCall = ''
  secondBestSigma = Inf
  secondBestClonality = 0
  for ( tryCall in allCalls()[allCalls() != call] ) {
    iscnv = isCNV(cluster, efs, callTofM(tryCall)['M'],
      callTofM(tryCall)['f'], callPrior(tryCall), bestSigma)
    if ( secondBestCall == '' || iscnv$sigma < secondBestSigma ) {
      secondBestCall = tryCall
      secondBestSigma = iscnv$sigma
      secondBestClonality = iscnv$clonality
    }
  }
  
  theoreticalError = abs(secondBestClonality - bestClonality)/
    (1 + sqrt(noneg(max(bestSigma, secondBestSigma)^2-bestSigma^2)))

  return(theoreticalError)
}


getBestSigma = function(cluster, call, efs) {
  bestSigma = Inf
  for ( tryCall in allCalls() ) {
    iscnv = isCNV(cluster, efs, callTofM(tryCall)['M'],
      callTofM(tryCall)['f'], callPrior(tryCall), 5)
    if ( iscnv$sigma < bestSigma ) {
      bestSigma = iscnv$sigma
    }
  }
  
  return(bestSigma)
}


makeTreeConsistent = function(cloneTree, cloneStories, storyList) {

  toRemove = findFirstUnitarityViolation(cloneTree, cloneStories, storyList)
  if ( length(toRemove) == 0 ) return(cloneTree)
  
  #remove story from tree.
  cloneStories = cloneStories[rownames(cloneStories) != toRemove,]
  storyList = storyList[names(storyList) != toRemove]

  cloneTree = findCloneTree(cloneStories)

  #recur until nothing needs to be removed
  cloneTree = makeTreeConsistent(cloneTree, cloneStories, storyList)
  
  return(cloneTree)
  
}

#looks recursively for unitarity violations in the clone tree, and returns the first main suspect of
#false clone that cause unitarity violation. returns empty vector if consistent with unitarity.
findFirstUnitarityViolation = function (cloneTree, cloneStories, storyList) {
  if (length(cloneTree) == 0) 
    return(c())
  for (name in names(cloneTree)) {
    if (length(cloneTree[[name]]) == 0) {
      next
    }
    whichViolating = whichUnitarityViolating(cloneStories[names(cloneTree[[name]]),], cloneStories[name,])
    if (length(whichViolating) > 0) {
      dodgyness = superFreq:::getDodgyness(storyList[whichViolating], 
        cloneStories[whichViolating,])
      suspect = whichViolating[which.max(dodgyness)]
      return(suspect)
    }
    deeperSuspect = findFirstUnitarityViolation(cloneTree[[name]], 
      cloneStories, storyList)
    if (length(deeperSuspect) > 0) {
      return(deeperSuspect)
    }
  }
  return(c())
}


#check which, if any, of the immediate subclones contribute to unitarity violation.
#returns a character vector of the contributing subclones, empty vector if consistent with unitarity.
whichUnitarityViolating = function(cloneStories, parentStory) {
  sampleTooHigh = superFreq:::colsums(superFreq:::noneg(cloneStories$stories - cloneStories$errors*1.5)) > (parentStory$stories+parentStory$errors*1.5)[1,]
  if ( !any(sampleTooHigh) ) return(c())
  contributing = superFreq:::rowsums(t(t(cloneStories$stories -
  cloneStories$errors*1.5 > 0)*sampleTooHigh)) > 0
  return(names(contributing)[contributing])
}


getDodgyness = function(storyList, stories) {
  dodgyness = sapply(names(storyList), function(name) {
    mutations = storyList[[name]]
    if ( name == 'germline' ) return(0)
    Nmut = length(mutations)
    Ncna = sum(grepl('^chr', mutations))
    Nindel = sum(grepl('[+-]', mutations))
    Nsnv = Nmut - Ncna - Nindel

    #few mutations overall
    ret = 1/Nmut
    
    #high indel/snv ratio
    ret = ret + Nindel/(Nsnv+Nindel+1)

    #close-by SNVs (less than 10kbp)
    if ( Nsnv > 1 ) {
      snvX = as.numeric(gsub('[A-Z]', '', mutations[!grepl('^chr', mutations) & !grepl('[+-]', mutations)]))
      snvX = sort(snvX)
      closeBy = snvX[2:length(snvX)] - snvX[1:(length(snvX)-1)] < 1e4
      ret = ret + noneg(sum(closeBy)/length(closeBy) - 0.1)
    }
    
    #high cna/snv ratio
    ret = ret + Ncna/(Ncna+Nsnv+1)

    #less dodgy if both CNA and SNV support
    ret = noneg(ret - (Ncna > 0 & Nsnv > 0)*0.2)

    #noise tends to be reasonably constant over samples
    highestLow = max(stories[name,]$stories - stories[name,]$errors, na.rm=T)
    lowestHigh = min(stories[name,]$stories + stories[name,]$errors, na.rm=T)
    significantChange = noneg(highestLow - lowestHigh)
    ret = ret + noneg(0.5 - significantChange)*2

    return(ret)
  })
  return(dodgyness)
}

clonesInTree = function(tree) {
  if ( length(tree) == 0 ) return('')

  namesOnThisLevel = names(tree)
  namesOnDeeperLevels = unlist(sapply(tree, clonesInTree))
  namesOnDeeperLevels = namesOnDeeperLevels[namesOnDeeperLevels != '']

  return(c(namesOnThisLevel, namesOnDeeperLevels))
}

renameClones = function(clusters) {
  before = clonesInTree(clusters$cloneTree)
  after = as.character(1:length(before))
  if ( before[1] == 'germline' )
    after = c('germline', 1:(length(before)-1))
  names(after) = before

  clusters$cloneTree = renameCloneTree(clusters$cloneTree, after)
  rownames(clusters$cloneStories) = after[rownames(clusters$cloneStories)]
  rownames(clusters$cloneStories$stories) = after[rownames(clusters$cloneStories$stories)]
  rownames(clusters$cloneStories$errors) = after[rownames(clusters$cloneStories$errors)]
  names(clusters$storyList) = after[names(clusters$storyList)]

  return(clusters)
}


renameCloneTree = function(tree, newNames) {
  if ( length(tree) == 0 ) return(tree)

  for ( i in 1:length(tree) ) {
    names(tree)[i] = newNames[names(tree)[i]]
    tree[[i]] = renameCloneTree(tree[[i]], newNames)
  }
  return(tree)
}
