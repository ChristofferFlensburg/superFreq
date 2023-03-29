

#' returns the 104 signature from a superFreq q variant object
#'
#' @import MutationalPatterns
#' @import BSgenome
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
get104profile = function(q, genome, somaticPcut = 0.5) {
  mx96 = get96signature(q=q, genome=genome, somaticPcut=somaticPcut)
  mxIndel = getIndelSignature(q, somaticPcut=somaticPcut)
  mx104 = rbind(mx96, mxIndel)
  return(mx104)
}


#get the 8 types of indel counts.
getIndelSignature = function(q, somaticPcut=0.5) {
  q = q[q$somaticP > somaticPcut,]
  ins1 = sum(grepl('^\\+', q$variant) & nchar(q$variant) == 2)
  ins2 = sum(grepl('^\\+', q$variant) & nchar(q$variant) %in% 3:4)
  ins4 = sum(grepl('^\\+', q$variant) & nchar(q$variant)  %in% 5:7)
  ins7p = sum(grepl('^\\+', q$variant) & nchar(q$variant) > 7)
  #the gsub throws a warning and returns NA when q$variant is not on '-[0-9]' format
  #those will turn FALSE from & grepl('^-', q$variant) though, so suppressing warnings here.
  del1 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) == 1 & grepl('^-', q$variant)))
  del2 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) %in% 2:3 & grepl('^-', q$variant)))
  del4 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) %in% 4:6 & grepl('^-', q$variant)))
  del7p = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) > 6 & grepl('^-', q$variant)))
  
  indelMx = matrix(c(ins1, ins2, ins4, ins7p, del1, del2, del4, del7p), ncol=1,
    dimnames=list(c('ins1', 'ins2', 'ins4', 'ins7p', 'del1', 'del2', 'del4', 'del7p'), 'name'))
  
  return(indelMx)
}

getRownamesFor96signatures = function() {
  bases_left = bases_right = c("A", "C", "G", "T")
  base_subs <- c("[C>A]", "[C>G]", "[C>T]", "[T>A]", "[T>C]", "[T>G]")
  full_context_poss <- vector("list", length(base_subs))
  for (i in seq_along(base_subs)) {
	  sub_context <- base_subs[[i]]
	  for (j in seq_len(1)) {
		  combi_tb <- tidyr::crossing(bases_left, sub_context, 
			  bases_right)
		  sub_context <- paste0(combi_tb$bases_left, combi_tb$sub_context, 
			  bases_right)
	  }
	  full_context_poss[[i]] <- sub_context
  }
  subs = unlist(full_context_poss)
  return(subs)
}


#get the 96 signature through mutationalPatterns
get96signature = function(q, genome, somaticPcut=0.5) {
  ref_genome = genomeToMPgenome(genome)
  qSom = q[q$somaticP > somaticPcut & q$variant %in% c('A', 'T', 'C', 'G'),]

  if ( nrow(qSom) == 0 ) {
  	ret = matrix(0,nrow=96, ncol=1)
  	colnames(ret) = 'temp'
  	rownames(ret) = getRownamesFor96signatures()
  	return(ret)
  }
  
  grl = qsToGRL(list('temp'=qSom), genome=genome)
  mx = mut_matrix(vcf_list = grl, ref_genome = ref_genome)  
  return(mx)
}

#converts a superFreq q data frame of variants to a GRanges made for MutationalPatterns
qToGR = function(q, genome) {
	#replace the legacy superFreq indel notation with what MutationalPAtterns expects
	q = superFreq:::updateDeletions(q, genome=genome)
	q = superFreq:::updateInsertions(q, genome=genome)
	chr = paste0('chr', xToChr(q$x, genome=genome))
	pos = xToPos(q$x, genome=genome)

	#QUAL from somaticP here is used with the usal phred conversion, pretending that somaticP is a probability
	gr = GRanges(seqnames=chr,
				 ranges=IRanges(start=pos, end=pos+nchar(q$reference)-1),
				 strand='*',
				 REF=q$reference,
				 ALT=q$variant,
				 QUAL = pmin(60, -log10(1 - q$somaticP) * 10),
				 FILTER = gsub("^$", "PASS", superFreq:::expandFlags(q$flag)),
				 INFO = paste0("DP=", q$cov, ";AF=", ifelse(q$cov > 0, q$var/q$cov, 0))
)	 
	
	chrL = chrLengths(genome)[gsub('chr', '', seqnames(seqinfo(gr)))]
	names(chrL) = seqnames(seqinfo(gr))
	seqlengths(seqinfo(gr)) = chrL
	
	GenomeInfoDb::genome(gr) = genome

	return(gr)
}
#converts a superFreq qs variants list to a GRanges list.
qsToGRL = function(qs, genome) {
	grl = lapply(qs, function(q) superFreq:::qToGR(q, genome=genome))
	return(grl)
	
}

#changes the legacy superFreq format for indels (such as "-3" for a 3bp deletion)
#to more up to date formats, where the reference entry has the deleted bases.
updateDeletions = function(q, genome) {
	dels = grep('-', q$variant)
	if ( length(dels) == 0 ) return(q)
	delSize = as.numeric(gsub('-', '', q$variant[dels]))
	delX = q$x[dels]
	ref_genome = superFreq:::genomeToMPgenome(genome)
	deletedSequence = as.character(getSeq(get(ref_genome), paste0('chr', xToChr(delX, genome)),
                         xToPos(delX, genome), xToPos(delX, genome)+delSize))

	q$x[dels] = q$x[dels]
	q$reference[dels] = deletedSequence
	q$variant[dels] = substr(deletedSequence, 1, 1)
    
    return(q)
}
updateInsertions = function(q, genome) {
	inss = grep('\\+', q$variant)
	if ( length(inss) == 0 ) return(q)
	q$variant[inss] = paste0(q$reference[inss], gsub('\\+', '', q$variant[inss]))
    
    return(q)
}


#translates the superFreq assembly names to the BSgenome assembly names
genomeToMPgenome = function(genome) {
  if ( genome == 'hg19' ) return('BSgenome.Hsapiens.UCSC.hg19')
  if ( genome == 'hg38' ) return('BSgenome.Hsapiens.UCSC.hg38')
  if ( genome == 'mm10' ) return('BSgenome.Mmusculus.UCSC.mm10')
  if ( genome == 'mm9' ) return('BSgenome.Mmusculus.UCSC.mm9')
  warning('Couldnt find genome')
  return('')
}

#plot the 104 profile
plot_104_profile = function (mut_matrix, ymax = 0, legend=T, main=paste0(colnames(mut_matrix), ' (', sum(mut_matrix), ' mutations)'), absolute=F, ...) {
  if ( any(class(mut_matrix) %in% c('numeric', 'integer')) & length(mut_matrix) == 104 ) mut_matrix = matrix(mut_matrix, ncol=1)
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/pmax(1,sum(x)))
  if ( absolute ) norm_mut_matrix = mut_matrix
  context = c(MutationalPatterns:::TRIPLETS_96, 'ins1', 'ins2', 'ins4', 'ins7', 'del1', 'del2', 'del4', 'del7')
  substitution = c(rep(MutationalPatterns:::SUBSTITUTIONS, each = 16), rep('indel',8))
  substring(context[1:96], 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = reshape2::melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if ( ymax==0 ) ymax = max(df3$value, 0.01)*1.05
  cols = superFreq:::mcri(c('red', 'green', 'magenta', 'blue', 'cyan', 'orange', 'violet'))
  a=barplot(df3$value, space = c(0.4,ifelse(df3$substitution[1:103] == df3$substitution[2:104], 0.4, 2)),
          col=cols[as.numeric(as.factor(df3$substitution))], ylim=c(-ymax*0.2,ymax), main=main, ...)

  baseToCol = c('C'=superFreq:::mcri('red'), 'G'=superFreq:::mcri('orange'), 'T'=superFreq:::mcri('blue'), 'A'=superFreq:::mcri('cyan'))
  base1 = gsub('\\..+$', '', df3$context)[1:96]
  base2 = gsub('>.+$', '', df3$substitution)[1:96]
  base3 = gsub('^.+]', '', df3$context)[1:96]
  xs = a[1:96,1]
  rect(xleft=xs-0.5, ybottom=-0.04*ymax, xright=xs+0.5, ytop=-0.01*ymax, border=F, col=baseToCol[base1])
  rect(xleft=xs-0.5, ybottom=-0.07*ymax, xright=xs+0.5, ytop=-0.04*ymax, border=F, col=baseToCol[base2])
  rect(xleft=xs-0.5, ybottom=-0.10*ymax, xright=xs+0.5, ytop=-0.07*ymax, border=F, col=baseToCol[base3])

  text(a[97:104,1], -ymax*0.05, df3$context[97:104], cex=0.6, srt=90)
  
  for ( subs in unique(df3$substitution) ) {
    x = mean(a[df3$substitution == subs,])
    text(x, -ymax*0.15, subs, cex=1, col=cols[as.numeric(as.factor(unique(df3$substitution)))][subs==unique(df3$substitution)])
  }

  if ( legend )
    legend('topleft', names(baseToCol), pch=15, col=baseToCol, pt.cex=2, cex=0.8)
    
  rownames(a) = context
  invisible(list(x=a, y=df3$value, ymax=ymax))
}


#plot the 104 profiles by sample
plot104profilesBySample = function(qs, plotDirectory, genome, somaticPcut=0.5) {
  patternDirectory = paste0(plotDirectory, '/mutationalPatterns')
  superFreq:::ensureDirectoryExists(patternDirectory)
  pdf(paste0(patternDirectory, '/bySample.pdf'), width=10, height=6)
  mxList = lapply(names(qs), function(sample) {
    mx104 = get104profile(q=qs[[sample]], genome=genome, somaticPcut=somaticPcut)
    superFreq:::plot_104_profile(mx104, main=paste0(sample, ' (', sum(mx104), ' mutations)'))
    return(mx104)
  })
  dev.off()
  mx = do.call(cbind, mxList)
  colnames(mx) = names(qs)
  invisible(mx)
}

#plot the 104 profiles by clone
plot104profilesByClone = function(storyList, qs, plotDirectory, genome) {
  patternDirectory = paste0(plotDirectory, '/mutationalPatterns')
  superFreq:::ensureDirectoryExists(patternDirectory)
  q = qs[[1]]
  pdf(paste0(patternDirectory, '/byClone.pdf'), width=10, height=6)
  mxList = lapply(names(storyList), function(clone) {
    mutations = storyList[[clone]]
    mutations = mutations[mutations %in% rownames(q)]
    if ( length(mutations) == 0 ) return()
    mx104 = superFreq:::get104profile(q=q[mutations,], genome=genome, somaticPcut=-1)
    colnames(mx104) = clone
    superFreq:::plot_104_profile(mx104, main=paste0('clone ', clone, ' (', sum(mx104), ' mutations)'))
    return(mx104)
  })
  dev.off()
  mx = do.call(cbind, mxList)
  invisible(mx)
}

#load the variants and stories from the R directory.
#calculate the profiles.
#plot the profiles by sample and clone.
#save the signatures to the Rdirectory.
loadAndPlotProfiles = function(Rdirectory, plotDirectory, genome) {
  storyFile = paste0(Rdirectory, '/stories.Rdata')
  variantFile = paste0(Rdirectory, '/allVariants.Rdata')
  if ( !all(file.exists(c(storyFile, variantFile))) ) {
    warning('Couldnt find data in ', Rdirectory, ' for mutational profiles. Could be wrong path or run hasnt finished.')
    return()
  }
  load(storyFile)
  stories = stories$stories
  load(variantFile)
  variants = allVariants$variants
  plotProfiles(plotDirectory=plotDirectory, stories=stories,
               variants=variants, genome=genome)
}

plotProfiles = function(Rdirectory, plotDirectory, stories, variants, genome, forceRedo=F) {
  catLog('Mutaional signatures of somatic mutations...')
  saveFile = paste0(Rdirectory, '/profiles.Rdata')
  if ( !file.exists(saveFile) | forceRedo ) {
    storyList = stories[[1]]$consistentClusters$storyList
    qs = variants$variants
    catLog('by sample..')
    sampleMx = plot104profilesBySample(qs=qs, plotDirectory=plotDirectory, genome=genome, somaticPcut=0.5)
    catLog('by clone..')
    cloneMx = plot104profilesByClone(storyList=storyList, qs=qs, plotDirectory=plotDirectory, genome=genome)
    profiles = list(sampleMx=sampleMx, cloneMx=cloneMx)
    catLog('saving results..')
    save(profiles, file=saveFile)
    catLog('done.\n')
  }
  else catLog('save file already there, moving along.\n')
}

#import functions for the preprocessed public signatures from COSMIC, TCGA and mutagens
importCOSMICsignatures = function(resourceDirectory='superFreqResources') {
  load(paste0(resourceDirectory, '/signatures/cosmicSignatures_v2_104.Rdata'))
  return(as.matrix(cosmicSignatures))
}
importTCGAsignatures = function(resourceDirectory='superFreqResources') {
  load(paste0(resourceDirectory, '/signatures/extendedPatternMx.Rdata'))
  return(extendedPatternMx)
}
#this one parses the public csv at https://data.mendeley.com/datasets/m7r4msjb4c/2, but the .Rdata is
#directly imported in superFreq, so this function is only an FYI for how we created the .Rdata
parseMutagenSignatures = function(txt='~/majewski/resources/public_signatures/Mutagen53_sub_signature.txt') {
  mutagenSignatures = read.table(txt, sep='\t', fill=T, stringsAsFactors=F)
  rownames(mutagenSignatures) = mutagenSignatures[,1]
  mutagenSignatures = mutagenSignatures[, -c(1,55)]
  mutagen = as.character(mutagenSignatures[1,])
  mutagenSignatures = mutagenSignatures[-1,]
  colnames(mutagenSignatures) = mutagen
  mutagenSignatures = rbind(mutagenSignatures, 'ins1'=0, 'ins2'=0, 'ins4'=0, 'ins7'=0, 'del1'=0, 'del2'=0, 'del4'=0, 'del7'=0)
  for ( i in 1:ncol(mutagenSignatures) ) mutagenSignatures[[i]] = as.numeric(mutagenSignatures[[i]])
  mutagenSignatures = as.matrix(mutagenSignatures)
  return(mutagenSignatures)
}
importMutagenSignatures = function(resourceDirectory='superFreqResources') {
  load(paste0(resourceDirectory, '/signatures/mutagenSignatures.Rdata'))
  return(as.matrix(mutagenSignatures))
}


plotSimilarSignatures = function(signature, resourceDirectory='superFreqResources', TCGA=T, COSMIC=T, mutagen=T, ignoreIndels=F, main='similar signatures', nLabel=20, nLabelTCGA=5, nLabelCOSMIC=5, nLabelMutagen=5, csvFile=NULL, ylab='mutation rate + 1', ...) {
  if ( !any(TCGA, COSMIC, mutagen) ) return()
  sources = c('TCGA', 'COSMIC', 'mutagen')[c(TCGA, COSMIC, mutagen)]
  sigList =
    list('TCGA'=superFreq:::importTCGAsignatures(resourceDirectory),
         'COSMIC'=superFreq:::importCOSMICsignatures(resourceDirectory),
         'mutagen'=superFreq:::importMutagenSignatures(resourceDirectory))[c(TCGA, COSMIC, mutagen)]
  sigMx = do.call(cbind, sigList)
  dataSource = unlist(lapply(sources, function(s) rep(s, ncol(sigList[[s]]))))

  if ( ignoreIndels ) {
    signature = c(signature[1:96], rep(0, 8))
    sigMx = rbind(sigMx[1:96,], matrix(0, ncol=ncol(sigMx), nrow=8))
  }
  
  sourceToCol = c('TCGA'=mcri('blue'), 'COSMIC'=mcri('red'), 'mutagen'=mcri('green'))
  sourceToLineCol = c('TCGA'=mcri('blue',al=0.5), 'COSMIC'=mcri('red',al=0.5), 'mutagen'=mcri('green',al=0.5))
  sourceToScatterCol = c('TCGA'='defaultBlue', 'COSMIC'='defaultRed', 'mutagen'='defaultGreen')
  sourceToCex = c('TCGA'=1.5, 'COSMIC'=2, 'mutagen'=2)
  
  sims = colSums(sigMx*signature)/sqrt(colSums(sigMx^2))/superFreq:::vecNorm(signature)
  sims[is.nan(sims)] = 0
  x = sims
  y = colSums(sigMx)+1 + superFreq:::noneg(1 + rnorm(ncol(sigMx), 0, .5))
  plot(x, y, type='n', xlab='cosine similarity', log='y', main=main, xlim=c(0,max(1, x*1.2)), ylab=ylab, ...)
  segments(1, 0.1, 1, 1e10, col='grey', lwd=0.5)
  for ( s in sources )
    plotColourScatter(x[dataSource==s], y[dataSource==s], cex=sourceToCex[s], col=sourceToScatterCol[s], add=T)
  
  label = rank(-sims) <= nLabel | (TCGA & rank(-sims-2*(dataSource=='TCGA')) <= nLabelTCGA) | (COSMIC & rank(-sims-2*(dataSource=='COSMIC')) <= nLabelCOSMIC) | (mutagen & rank(-sims-2*(dataSource=='mutagen')) <= nLabelMutagen)
  labelPos = spreadPositions2D(x, log10(y), 0.07, 0.2, label, verbose=F)   #add to superFreq package, or external dependency
  labelText = gsub('TCGA\\.TCGA-', '', gsub('COSMIC\\.Signature\\.', 'COSMIC-', gsub('mutagen\\.X\\.', '', colnames(sigMx))))
  text(labelPos[1,], 10^labelPos[2,], colnames(sigMx)[label], col=sourceToCol[dataSource[label]], cex=0.7, font=2)
  diff = rbind(x[label], log10(y[label])) - labelPos
  dist = sqrt((labelPos[1,]-x[label])^2+(labelPos[2,]-y[label])^2)
  segments(labelPos[1,] + 0.1*diff[1,], 10^(labelPos[2,] + 0.1*diff[2,]), x[label], y[label], col=sourceToLineCol[dataSource[label]], lwd=0.5)

  if ( !is.null(csvFile) ) {
    dat = data.frame(participant=colnames(sigMx), similarity=sims, mutation_rate=colSums(sigMx))
    dat = dat[order(-dat$similarity),]
    write.table(dat, sep=',', file=csvFile, quote=F,row.names=F)
  }
  
  invisible(list('similarity'=sims, 'rates'=colSums(extendedPatternMx)))
}


#Takes a set of points with coordinates x, and y, of which the logical vector toLabel
#denotes the points to be labeled. xR and yR are the size (the radius) of the labels,
#used for collision ellipse.
spreadPositions2D = function(x, y, xR, yR, toLabel, verbose=T) {
  #starting condition
  labelPos = sapply(which(toLabel), function(i) c('x'=x[i], 'y'=y[i]))
  origin = labelPos
  allX = c(labelPos[1,], x)
  allY = c(labelPos[2,], y)
  relDist = apply(labelPos, 2, function(pos) (pos[1] - allX)^2/xR^2 + (pos[2] - allY)^2/yR^2)
  overlap = relDist < 1  #maybe should use square collision box instead?
  overlap[row(overlap) == col(overlap)] = 0

  #randomly move away from locally densest area until not overlapping 
  while ( sum(overlap) > 0 ) {
    labelPos = sapply(1:ncol(labelPos), function(id) {
      dist = relDist[,id]
      dist[id] = 100
      closest = (dist < runif(1,1,6))
      posOri = labelPos[,id]
      if ( min(dist) >= 1 ) return(posOri)
      pos0 = c('x'=mean(allX[closest]), 'y'=mean(allY[closest]))
      diff = posOri - pos0
      if ( diff[1] == 0 & diff[2] == 0 ) diff = c('x'=0, 'y'=1)
      newPos = posOri + runif(1,0.1,0.5)*diff/sqrt((diff[1]/xR)^2 + (diff[2]/yR)^2)
      return(newPos)
    })

    allX = c(labelPos[1,], x)
    allY = c(labelPos[2,], y)
    relDist = apply(labelPos, 2, function(pos) (pos[1] - allX)^2/xR^2 + (pos[2] - allY)^2/yR^2)
    overlap = relDist < 1
    overlap[row(overlap) == col(overlap)] = 0

    if ( verbose ) cat(colSums(overlap), '\n')
  }

  #tighten rubber band towards starting point
  tight = F
  while ( !tight ) {
    tight=T
    for ( id in 1:ncol(labelPos) ) {
      diff = origin[,id] - labelPos[,id]
      dist = sqrt((labelPos[1,]-x[toLabel])^2+(labelPos[2,]-y[toLabel])^2)
      newPos = labelPos[,id] + 0.1*diff
      
      dists = sqrt((newPos[1] - allX)^2/xR^2 + (newPos[2] - allY)^2/yR^2)
      dists[id] = 100
      if ( min(dists) >= 0.99 ) {
        if ( verbose ) cat('tightening', id, '\n')
        tight=F
        labelPos[,id] = newPos
        allX = c(labelPos[1,], x)
        allY = c(labelPos[2,], y)
      }
    }
  }

  return(labelPos)
}




