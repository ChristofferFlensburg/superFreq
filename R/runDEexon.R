
#This function take a set of bamfiles, a set of normal bamfiles, and capture regions as input.
#The function runs differential coverage one each sample vs the pool of normals using limma-voom.
#the counts are loess and binding strength corrected.

#' Run differential coverage analysis.
#'
#' @param bamFiles vector of paths to the bamfiles to be analysed
#' @param names vector of names of the samples associated with the bam files.
#' @param externalNormalBams vector of paths to the pool of normal samples.
#' @param captureRegions GRanges capture regions, as you get from importCaptureRegions.
#' @param Rdirectory The directory to save or load the output to/from.
#' @param plotDirectory The directory to plot results and diagnosis to.
#' @param normalRdirectory The directory to save or load the pool of normal counts to/from.
#' @param settings Settings as passed to analyse(...).
#' @param genome Character string of the genome, such as 'hg19', 'hg38' or 'mm10'. Defaults to 'hg19'.
#' @param cpus The number of parallel cpus to be used. Defaults to 1.
#' @param forceRedoFit Boolean, forcing the differential analysis to be redone even if saved data is
#'                     present. Defaults to FALSE.
#' @param forceRedoCount Boolean, forcing the counting to be redone even if saved data is
#'                     present. Defaults to FALSE.
#' @param forceRedoNormalCount Boolean, forcing the counting of normal samples to be redone even
#'                     if saved data is present. Defaults to FALSE.
#'
#' @details Reads in the bed file into GRanges format, with extra columns for GC and binding strength.
#'
#' @importFrom BiocGenerics strand start end
#' @importFrom Rsubread featureCounts
#' @importFrom limma makeContrasts voom lmFit eBayes contrasts.fit
#' @importFrom parallel mclapply
#'
runDE = function(bamFiles, names, externalNormalBams, captureRegions, Rdirectory, plotDirectory, normalRdirectory,
  settings=list(), genome='hg19', cpus=1,
  forceRedoFit=F, forceRedoCount=F, forceRedoNormalCount=F) {
  catLog('Starting differential coverage analysis by sample.\n')

  fitSaveFile = paste0(Rdirectory, '/fit.Rdata')
  if ( file.exists(fitSaveFile) & !forceRedoFit & !forceRedoCount & !forceRedoNormalCount ) {
    catLog('Loading saved differential coverage results.\n')
    load(file=fitSaveFile)
    catLog('Loaded saved differential coverage results\n')
    return(fit)
  }

  catLog('Preparing capture regions for featureCounts..')
  captureAnnotation = try(captureRegionToAnnotation(captureRegions))
  if ( class('captureAnnotation') == 'try-error' ) stop('Error in captureRegionToAnnotation.')
  if ( !('GeneID' %in% colnames(captureAnnotation)) ) stop('captureAnnotation does not have a GeneID.')
  catLog('done.\n')
  
  fCsSaveFile = paste0(Rdirectory, '/fCsExon.Rdata')
  if ( !file.exists(fCsSaveFile) | forceRedoCount ) {
    catLog('Counting reads over capture regions.\n')
    fCsExon = try(featureCounts(bamFiles, annot.ext=captureAnnotation, useMetaFeatures=F,
      allowMultiOverlap=T, isPairedEnd=T, minMQS=10, nthreads=cpus))
    if ( class(fCsExon) != 'list' ) {
      catLog('Error in featureCounts.\nInput was\nbamFiles:', bamFiles,
             '\ncaptureAnnotation[1:10,]:', as.matrix(captureAnnotation[1:10,]), '\n')
      stop('Error in featureCounts.')
    }
    catLog('Got a sample count matrix of size', dim(fCsExon$counts), ', with total counts:\n')
    for ( col in colnames(fCsExon$counts) ) {
      catLog(col, ': ', sum(fCsExon$counts[,col]), '\n')
    }
    colnames(fCsExon$counts) = names
    catLog('Saving sample counts to ', fCsSaveFile, '..', sep='')
    save(fCsExon, file=fCsSaveFile)
    catLog('done.\n')
  }
  else {
    catLog('Loading counts from file..')
    load(fCsSaveFile)
    catLog('done.\n')
    catLog('Loaded counts of dimension', dim(fCsExon$counts), '\n')
  }

  normalFCsSaveFile = paste0(normalRdirectory, '/normalFCsExon.Rdata')
  if ( !file.exists(normalFCsSaveFile) | forceRedoNormalCount ) {
    catLog('Counting normal reads over capture regions.\n')
    normalFCsExon = try(featureCounts(externalNormalBams, annot.ext=captureAnnotation, useMetaFeatures=F,
      allowMultiOverlap=T, isPairedEnd=T, minMQS=10, nthreads=cpus))
    if ( class(normalFCsExon) != 'list' ) {
      catLog('Error in featureCounts.\nInput was\nexternalNormalBams:', externalNormalBams,
             '\ncaptureAnnotation[1:10,]:', as.matrix(captureAnnotation[1:10,]), '\n')
      stop('Error in featureCounts of normals.')
    }
    catLog('Got a normals count matrix of size', dim(normalFCsExon$counts), ', with total counts:\n')
    for ( col in colnames(normalFCsExon$counts) ) {
      catLog(col, ': ', sum(normalFCsExon$counts[,col]), '\n')
    }
    colnames(normalFCsExon$counts) = names(externalNormalBams)
    catLog('Saving normals counts to ', normalFCsSaveFile, '..', sep='')
    save(normalFCsExon, file=normalFCsSaveFile)
    catLog('done.\n')
  }
  else {
    catLog('Loading normals counts from file..')
    load(normalFCsSaveFile)
    catLog('done.\n')
    catLog('Loaded normal counts of dimension', dim(normalFCsExon$counts), '\n')
    if ( length(externalNormalBams) != ncol(normalFCsExon$counts) ) {
      warning('Saved normal counts dont match the bam files curently in the normal directory. Please rerun with \"forceRedo$forceRedoNormalCount = T\" to redo the normal counts on the current bam files.')
    }
  }
  catLog('Merging sample and normals counts..')
  if ( nrow(fCsExon$counts) != nrow(normalFCsExon$counts) ) {
    warning('Seems like the normals and the analysed samples are not counted over the same capture regions,\n',
            'or at least they dont have the same number of counts.\n',
           'Try redoing the counts for both by setting forceRedo$forceRedoCount and forceRedo$forceRedoNormalCount to TRUE.')
  }
  counts = cbind(fCsExon$counts, normalFCsExon$counts)
  annotation = fCsExon$annotation
  catLog('done.\n')

  if ( any(colSums(counts)==0) ) {
    warning('0 reads counted in sample(s): ', colnames(counts)[colSums(counts)==0], '. This will cause nonsensical CNA calls. These counts are from featureCounts.')
    counts[,colSums(counts)==0] = 1
  }

  catLog('Determining sex..')
  xes = grep('^X', annotation$Chr)
  yes = grep('^Y', annotation$Chr)
  if ( length(xes)==0 | length(yes)==0 ) {
    catLog('Capture regions not present in both chromosome X and Y, skip sex matching.\n')
    sex = rep('female', ncol(counts))
    names(sex) = colnames(counts)
  }
  else {
    maleScore =
      colsums(10*libNorm(counts)[yes,])/sum(annotation$Length[yes]+300) -
        colsums(libNorm(counts)[xes,])/sum(annotation$Length[xes]+300)
    sex = ifelse(maleScore > 0, 'male', 'female')
    names(sex) = colnames(counts)
    catLog('done.\nSAMPLE', 'SCORE', 'SEX', '\n', sep='   ')
    for ( i in 1:length(sex) ) catLog(names(sex)[i], maleScore[i], sex[i], '\n', sep='   ')
  }

  catLog('Setting up design matrix for linear analysis..')
  group = c(names, paste0(rep('normal', length(externalNormalBams))))
  design = model.matrix(~0+group)
  colnames(design) = gsub('^group', '', gsub('^sex', '',colnames(design)))
  contrastList = c(lapply(names, function(name) paste0(name, '-normal')), list(levels=colnames(design)))
  contrasts = do.call(makeContrasts, contrastList)
  catLog('done.\n')
  catLog('Design matrix is\n', colnames(design), '\n', sep='   ')
  for ( row in rownames(design) ) catLog(design[row,], '\n', sep=' ')

  #MA plots
  diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
  if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
  catLog('Making MA plots of coverage in diagnostics directory..')
  normalMean = rowMeans(counts[,group=='normal'])
  MAdirectory = paste0(diagnosticPlotsDirectory, '/MAplots/')
  if ( !file.exists(MAdirectory) ) dir.create(MAdirectory)
  for ( col in colnames(counts) ) {
    catLog(col, '..', sep='')
    plotfile = paste0(MAdirectory, col, '.png')
    if ( !file.exists(plotfile) | forceRedoFit ) {
      png(plotfile, height=10, width=20, res=144, unit='in')
      plotMA(counts[,col], normalMean, loess=T, span=0.5, verbose=F, main = paste0(col, ' vs normals (before loess normalisation)'))
      dev.off()
    }
  }
  catLog('done.\n')

  countsOri = counts

  if ( getSettings(settings, 'sexCorrection') ) {
    catLog('Correcting for sex chromsome coverage..')
    sexCorrectionFactor = list()
    #multiply male coverage by 2 over the x chromsome in the normal samples
    sexRows = 1:nrow(counts) %in% c(xes, yes)
    for (i in 1:ncol(counts)) {
      col = colnames(counts)[i]
      if ( sex[col] == 'male' ) {
        nonSexReads = sum(counts[!sexRows,col])
        sexReads = sum(counts[sexRows,col])
        normalisationFactor = (nonSexReads+sexReads)/(nonSexReads+2*sexReads)
        sexCorrectionFactor[[i]] = ifelse(sexRows, 2*normalisationFactor, normalisationFactor)
      }
      else sexCorrectionFactor[[i]] = rep(1, nrow(counts))
    }
    sexCorrectionFactor = do.call(cbind, sexCorrectionFactor)
    countsSex = counts*sexCorrectionFactor
    catLog('done.\n')
  }
  else {
    catLog('Skipping sex correction.\n')
    countsSex = counts
  }


  #loess normalise normals, then match the non-normal M-A dependence to the normals.
  #This is not needed for most samples, but will really improve accuracy when there is an M-A bias.
  if ( getSettings(settings, 'MAcorrection') ) {
    catLog('Loess normalising counts to normals..')
    counts = countsSex
    counts[,group=='normal'] = loessNormAll(counts[,group=='normal',drop=F], span=0.5)
    counts[,group=='normal'] = loessNormAll(counts[,group=='normal',drop=F], span=0.5)
    counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal',drop=F], counts[,group=='normal',drop=F], span=0.5)
    counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal',drop=F], counts[,group=='normal',drop=F], span=0.5)
    loessCorrectionFactor = (0.5+counts)/(0.5+countsSex)
    countsSexLo = counts
    catLog('done.\n')
  }
  else {
    catLog('Skipping loess normalising.\n')
    counts = countsSex
    countsSexLo = counts
  }

  #local binding strength correction
  if ( getSettings(settings, 'GCcorrection') ) {
    catLog('Correcting for binding strength bias..')
    BSdirectory = paste0(diagnosticPlotsDirectory, '/BSplots/')
    if ( !file.exists(BSdirectory) ) dir.create(BSdirectory)
    dn = as.numeric(captureRegions$dn)
    counts = countsSexLo
    w = pmin(2, pmax(0.1, sqrt(rowMeans(countsOri)/100)))
    maxCorrection = 10
    BScorrectionFactor =
      lapply(1:ncol(counts), function(i) {
        col = colnames(counts)[i]
        catLog(col, '..', sep='')
        LFC = log((1+counts[,col])/(annotation$Length))
        lo = loess(cov~dn, data=data.frame(cov=LFC, dn=dn),
          weights=w, span=0.3, family='symmetric', degrees=1, trace.hat='approximate')
        los = predict(lo, dn)
        biasFactor = exp(los)
        correctionFactor = (1/biasFactor)/mean(1/biasFactor)
        correctionFactor = correctionFactor*sum(counts[,col])/sum(counts[,col]*correctionFactor)
        #find normalisation such that no count is corrected more than a factor 2 and maintaining normalisation
        trialN = (500:2000)/1000
        trialTotCount = sapply(trialN, function(N)
          sum(counts[,col]*pmax(1/maxCorrection, pmin(maxCorrection,N*correctionFactor))))
        totCount = sum(counts[,col])
        diff = trialTotCount - totCount
        N = trialN[which(abs(diff) == min(abs(diff)))[1]]
        correctionFactor = pmax(1/maxCorrection, pmin(maxCorrection,N*correctionFactor))
        plotfile = paste0(BSdirectory, col, '.png')
        if ( !file.exists(plotfile) | forceRedoFit ) {
          use = dn > 0
          ylim = quantile(LFC[use], probs=c(0.01, 0.99))
          png(plotfile, height=10, width=20, res=144, unit='in')
          plotColourScatter(dn[use], LFC[use], cex=w[use], ylim=ylim,
                            main=col, xlab='binding strength', ylab='log(reads/bp)')
        lines((500:1500)/100, predict(lo, (500:1500)/100), lwd=5, col=mcri('orange'))
          legend('topleft', 'loess fit', lwd=5, col=mcri('orange'), bg='white')
          dev.off()
        }
        plotfile = paste0(BSdirectory, col, '.overGenome.png')
        if ( !file.exists(plotfile) | forceRedoFit ) {
          ylim = quantile(LFC, probs=c(0.01, 0.99))+c(0.5,0.5)
          png(plotfile, height=10, width=20, res=144, unit='in')
          plotColourScatter(annotationToX(annotation, genome), LFC, cex=w, ylim=ylim, xlab='genome',
                            ylab='~log(1+read depth)')
          points(annotationToX(annotation, genome), (dn-mean(dn))/1.5 + mean(ylim), cex=w/2, pch=16,
                 col=mcri('orange', 0.5))
          addChromosomeLines(ylim=ylim, col=mcri('green'), genome=genome)
          legend('bottomright', c('binding strength'), pch=16, col=mcri('orange'), bg='white')
          dev.off()
        }
        plotfile = paste0(BSdirectory, col, '.overGenome.BScorrected.png')
        if ( !file.exists(plotfile) | forceRedoFit ) {
          ylim = quantile(LFC + log(correctionFactor), probs=c(0.01, 0.99))+c(0.5,0.5)
          png(plotfile, height=10, width=20, res=144, unit='in')
          plotColourScatter(annotationToX(annotation, genome), LFC, cex=w, ylim=ylim, xlab='genome',
                            ylab='~log(1+read depth)')
          points(annotationToX(annotation, genome), (dn-mean(dn))/1.5 + mean(ylim), cex=w/2, pch=16,
                 col=mcri('orange', 0.5))
          addChromosomeLines(ylim=ylim, col=mcri('green'), genome=genome)
          legend('bottomright', c('binding strength'), pch=16, col=mcri('orange'), bg='white')
          dev.off()
        }
        return(correctionFactor)
      })
    BScorrectionFactor = do.call(cbind, BScorrectionFactor)
    countsSexLoBS = countsSexLo*BScorrectionFactor
    #maintain total count for each gene. Removed, as it increases interaction between GC and loess corrections.
    #countsSexLoBS = countsSexLoBS*(1+rowSums(countsOri))/(1+rowSums(countsSexLoBS))
    catLog('done!\n')
  }
  else {
    catLog('Skipping correction for GC bias.\n')
    countsSexLoBS = countsSexLo
  }

  

  #binding strength correcting with running average
  if ( getSettings(settings, 'GCregionCorrection') ) {
    catLog('Correcting for regional binding strength bias..')
    x = annotationToX(annotation, genome)
    counts = countsSexLoBS
    w = pmin(2, pmax(0.1, sqrt(rowMeans(countsOri)/100)))
    regionIs = mclapply(x, function(X) which(abs(x-X) < 5e5), mc.cores=cpus)
    regionDN = sapply(regionIs, function(is) sum(dn[is]*w[is])/sum(w[is]))
    maxCorrection = 10
    DNregionCorrectionFactor =
      mclapply(1:ncol(counts), function(i) {
        col = colnames(counts)[i]
        catLog(col, '..', sep='')
        LFC = log((1+counts[,col])/(annotation$Length))
        smoothLFC = sapply(regionIs, function(is) sum(LFC[is]*w[is])/sum(w[is]))
        lo = loess(LFC~regionDN, data=data.frame(LFC=smoothLFC, regionDN=regionDN),
          weights=w, span=0.3, degrees=1, trace.hat='approximate', family='symmetric')
        los = predict(lo, regionDN)
        biasFactor = exp(los)
        correctionFactor = (1/biasFactor)/mean(1/biasFactor)
        correctionFactor = correctionFactor*sum(counts[,col])/sum(counts[,col]*correctionFactor)
        #find normalisation such that no count is corrected more than a factor 10 and maintaining normalisation
        trialN = (500:2000)/1000
        trialTotCount = sapply(trialN, function(N)
          sum(counts[,col]*pmax(1/maxCorrection, pmin(maxCorrection,N*correctionFactor))))
        totCount = sum(counts[,col])
        diff = trialTotCount - totCount
        N = trialN[which(abs(diff) == min(abs(diff)))[1]]
        correctionFactor = pmax(1/maxCorrection, pmin(maxCorrection,N*correctionFactor))
        plotfile = paste0(BSdirectory, col, '.region.png')
        if ( !file.exists(plotfile) | forceRedoFit ) {
          ylim = quantile(smoothLFC, probs=c(0.01, 0.99))
          png(plotfile, height=10, width=20, res=144, unit='in')
          plotColourScatter(regionDN, smoothLFC, cex=w, main=col, ylim=ylim,
                            xlab='region binding strength', ylab='region mean log(reads/bp)')
          lines((500:1500)/100, predict(lo, (500:1500)/100), lwd=5, col=mcri('orange'))
          legend('topleft', 'loess fit', lwd=5, col=mcri('orange'))
          dev.off()
        }
        plotfile = paste0(BSdirectory, col, '.region.overGenome.png')
        if ( !file.exists(plotfile) | forceRedoFit ) {
          ylim = quantile(LFC, probs=c(0.01, 0.99))+c(0.5,0.5)
          png(plotfile, height=10, width=20, res=144, unit='in')
          plotColourScatter(annotationToX(annotation, genome), LFC, cex=w, ylim=ylim,
                            xlab='genome', ylab='~log(1+read depth)')
          points(annotationToX(annotation, genome), (regionDN-mean(regionDN))/1.5 + mean(ylim), cex=w/2,
                 pch=16, col=mcri('orange', 0.5))
          addChromosomeLines(ylim=ylim, col=mcri('green'), genome=genome)
          legend('bottomright', c('binding strength'), pch=16, col=mcri('orange'), bg='white')
          dev.off()
        }
        return(correctionFactor)
      }, mc.cores=cpus)
    DNregionCorrectionFactor = do.call(cbind, DNregionCorrectionFactor)
    countsSexLoBSregion = countsSexLoBS*DNregionCorrectionFactor
    #maintain total count for each gene. Removed, as it increases interaction between GC and loess corrections.
    #countsSexLoBSregion = countsSexLoBSregion*(1+rowSums(countsOri))/(1+rowSums(countsSexLoBSregion))
    catLog('done!\n')
  }
  else {
    catLog('Skipping correction for regional GC bias.\n')
    countsSexLoBSregion = countsSexLoBS
  }


  #loess correct the BS corrected counts
  if ( getSettings(settings, 'MAcorrection') ) {
    catLog('Second round of loess normalisation..')
    counts = countsSexLoBSregion
    counts[,group=='normal'] = loessNormAll(counts[,group=='normal',drop=F], span=0.5)
    counts[,group=='normal'] = loessNormAll(counts[,group=='normal',drop=F], span=0.5)
    counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal',drop=F], counts[,group=='normal',drop=F], span=0.5)
    counts[,group!='normal'] = loessNormAllToReference(counts[,group!='normal',drop=F], counts[,group=='normal',drop=F], span=0.5)
    loessCorrectionFactor2 = (0.5+counts)/(0.5+countsSexLoBS)
    countsSexLoBSLo = counts
    catLog('done.\n')
  }
  else {
    catLog('Skipping second round of loess normalisation.\n')
    counts = countsSexLoBSregion
    countsSexLoBSLo = counts
  }



  if ( getSettings(settings, 'sexCorrection') ) {
    catLog('Returning sex effects to non-normal samples..')
    #divide back male coverage by 2 over the x chromsome in the (non-normal) samples
    returnSampleSexCorrectionFactor = list()
    sexRows = 1:nrow(counts) %in% c(xes, yes)
    for (i in 1:ncol(counts)) {
      col = colnames(counts)[i]
      if ( sex[col] == 'male' & i <= length(names) ) {
        returnSampleSexCorrectionFactor[[i]] = ifelse(sexRows, 1/2, 1)
      }
      else returnSampleSexCorrectionFactor[[i]] = rep(1, nrow(counts))
    }
    returnSampleSexCorrectionFactor = do.call(cbind, returnSampleSexCorrectionFactor)
    countsSexLoBSLoSex = countsSexLoBSLo*returnSampleSexCorrectionFactor
    counts = countsSexLoBSLoSex
    catLog('done.\n')
  }
  else {
    catLog('Skipping second round of sex correction.\n')
    sexRows = 1:nrow(counts) %in% c(xes, yes)
    countsSexLoGCLoSex = countsSexLoGCLo
    counts = countsSexLoGCLoSex
  }


  improvement = c()
  breakImprovement = c()
  noise = c()
  breaks = c()
  normalMean = rowMeans(counts[,group=='normal'])
  normalMeanOri = rowMeans(countsOri[,group=='normal'])
  for ( col in 1:ncol(counts) ) {
    first = 1:(nrow(counts)-1)
    second = 2:nrow(counts)
    LFC = log((1+counts[,col])/(1+normalMean))
    w = 1/sqrt(1+counts[,col]) + 1/sqrt(1+normalMean)/sqrt(sum(group=='normal'))
    W = sqrt(w[first]^2+w[second]^2)
    change = abs(LFC[second]-LFC[first])/W
    LFCori = log((1+countsOri[,col])/(1+normalMeanOri))
    wOri = 1/sqrt(1+countsOri[,col]) + 1/sqrt(1+normalMeanOri)/sqrt(sum(group=='normal'))
    Wori = sqrt(wOri[first]^2+wOri[second]^2)
    changeOri = abs(LFCori[second]-LFCori[first])/Wori
    noise = c(noise, sum(changeOri > 1)/length(change))
    breaks = c(breaks, sum(changeOri > 10)/length(change))
    improvement = c(improvement, sum(change > 1)/sum(changeOri > 1))
    breakImprovement = c(breakImprovement, sum(change > 10)/sum(changeOri > 10))
  }
  names(improvement) = names(noise) = names(breakImprovement) = names(breaks) = colnames(counts)
  catLog('Corrections changes average LFC between sample by a factor ', round(mean(improvement) ,4),'.\n')

  #MA plots after loess and BS correction
  catLog('Making MA plots of coverage after loess and BS correction in diagnostics directory..')
  normalMean = rowMeans(counts[,group=='normal'])
  MAdirectory = paste0(diagnosticPlotsDirectory, '/MAplotsAfter/')
  if ( !file.exists(MAdirectory) ) dir.create(MAdirectory)
  for ( col in colnames(counts) ) {
    catLog(col, '..', sep='')
    plotfile = paste0(MAdirectory, col, '.exon.png')
    if ( !file.exists(plotfile) | forceRedoFit ) {
      png(plotfile, height=10, width=20, res=144, unit='in')
      plotMA(counts[,col], normalMean, loess=T, span=0.5, verbose=F, main = paste0(col, ' vs normals (after corrections)'))
      dev.off()
    }
  }
  catLog('done.\n')  


  #remove Y counts if no normal sample is male.
  if ( !any(sex[group=='normal'] == 'male') & length(yes) > 0 ) {
    counts = counts[-yes,]
    countsOri = countsOri[-yes,]
    annotation = annotation[-yes,]
    xes = grep('X$', annotation$Chr)
    yes = grep('Y$', annotation$Chr)
  }

  #merge exon counts to gene counts
  geneCountsOri = aggregate(x = countsOri, list(rownames(countsOri)), FUN = 'sum')
  rownames(geneCountsOri) = geneCountsOri[,1]
  geneCountsOri = geneCountsOri[,-1]
  geneCounts = aggregate(x = counts, list(rownames(counts)), FUN = 'sum')
  rownames(geneCounts) = geneCounts[,1]
  geneCounts = geneCounts[,-1]
  
  x = annotationToX(annotation, genome=genome)
  geneX = aggregate(x, list(rownames(counts)), FUN='min')[,2]
  geneOrder = order(geneX)
  geneX = geneX[geneOrder]
  geneCountsOri = geneCountsOri[geneOrder,]
  geneCounts = geneCounts[geneOrder,]
  geneXes = which(aggregate(grepl('^X', annotation$Chr), list(rownames(counts)), FUN = 'any')[geneOrder,2])
  geneYes = which(aggregate(grepl('^Y', annotation$Chr), list(rownames(counts)), FUN = 'any')[geneOrder,2])

  catLog('Running voom on exons..')
  png(paste0(diagnosticPlotsDirectory, '/voomVariance.exon.png'), height=10, width=20, res=144, unit='in')
  voomWeights = voomWithQualityWeights(countsOri, design=design, plot=T)
  voomWeights$weights[yes,design[,'normal']==1 & sex == 'female'] = 0
  dev.off()
  voomCounts = voomWithQualityWeights(counts, design=design, plot=F)
  voomCombined = voomCounts
  voomCombined$weights = voomWeights$weights
  catLog('limma..')
  exonFit = lmFit(voomCombined, design=design)
  exonFit = eBayes(exonFit)
  exonFit = contrasts.fit(exonFit, contrasts)
  exonFit = eBayes(exonFit)
  catLog('XRank..')
  exonFit = XRank(exonFit, plot=F, cpus=cpus, verbose=T, keepPosterior=F)
  exonFit$x = annotationToX(annotation, genome=genome)
  x1x2 = annotationToX1X2(annotation, genome=genome)
  exonFit$x1 = x1x2[,1]
  exonFit$x2 = x1x2[,2]
  exonFit$chr = annotationToChr(annotation)
  exonFit$longNames = annotation$GeneID
  catLog('done.\n')

  #remember stats from counting
  catLog('Importing stats about feature counts..')
  tots = colSums(fCsExon$stat[,-1,drop=F])
  assigned = fCsExon$stat[1,-1]
  names(tots) = names(assigned) = names
  exonFit$totNReads = tots
  exonFit$assignedNReads = assigned
  fitSex = sex[1:ncol(exonFit)]
  names(fitSex) = colnames(exonFit)
  exonFit$sex = fitSex
  catLog('done.\n')

  catLog('Running voom on genes..')
  png(paste0(diagnosticPlotsDirectory, '/voomVariance.exon.png'), height=10, width=20, res=144, unit='in')
  voomWeights = voomWithQualityWeights(geneCountsOri, design=design, plot=T)
  voomWeights$weights[geneYes,design[,'normal']==1 & sex == 'female'] = 0
  dev.off()
  voomCounts = voomWithQualityWeights(geneCounts, design=design, plot=F)
  voomCombined = voomCounts
  voomCombined$weights = voomWeights$weights
  catLog('limma..')
  fit = lmFit(voomCombined, design=design)
  fit = eBayes(fit)
  fit = contrasts.fit(fit, contrasts)
  fit = eBayes(fit)
  catLog('XRank..')
  fit = XRank(fit, plot=F, cpus=cpus, verbose=T, keepPosterior=F)
  fit$x1 = aggregate(exonFit$x1, list(rownames(counts)), FUN='min')[geneOrder,2]
  fit$x2 = aggregate(exonFit$x2, list(rownames(counts)), FUN='max')[geneOrder,2]
  fit$x = geneX
  fit$chr = xToChr(fit$x, genome)
  fit$sex = fitSex
  fit$longNames = rownames(geneCounts)
  catLog('done.\n')

  fit = list('fit'=fit, 'exonFit'=exonFit)
  
  catLog('Saving fit..')
  save(fit, file=fitSaveFile)
  catLog('done.\n')

  catLog('Returning fit of dimension', dim(fit$fit), 'and', dim(fit$exonFit),'\n')
  return(fit)
}










#helper function transforming GRanges capture regions to a format that can be inserted into featureCounts.
captureRegionToAnnotation = function(cR) {
  ret = data.frame(GeneID = as.character(cR$region), Start = start(cR), End = end(cR), Strand=strand(cR), Chr = seqnames(cR), stringsAsFactors=F)
  return(ret)
}

#helper functions doing loess normalisation from MA plots.
loessNormAll = function(counts, ...) {
  newCounts = do.call(cbind, lapply(1:ncol(counts), function(col) loessNorm(counts[,col], counts[,-col], ...)[,1]))
  colnames(newCounts) = colnames(counts)
  rownames(newCounts) = rownames(counts)
  return(newCounts)
}
loessNormAllToReference = function(counts, reference, ...) {
  newCounts = do.call(cbind, lapply(1:ncol(counts), function(col) loessNormToReference(counts[,col], reference, ...)))
  colnames(newCounts) = colnames(counts)
  rownames(newCounts) = rownames(counts)
  return(newCounts)
}
loessNorm = function(counts1, counts2, span=0.5) {
  if ( is.matrix(counts1) ) x = rowSums(1+counts1)
  else x = 1 + counts1
  if ( is.matrix(counts2) ) y = rowSums(1+counts2)
  else y = 1 + counts2
  M = log(x/y)
  A = log(x*y)/2
  lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
  loM = predict(lo, A)
  loM = loM - sum(loM*exp(A))/sum(exp(A))

  counts1 = counts1*exp(-loM/2)
  counts2 = counts2*exp(loM/2)
  return(cbind(counts1, counts2))
}
loessNormToReference = function(counts1, reference, span=0.5) {
  if ( is.matrix(counts1) ) x = rowSums(1+counts1)
  else x = 1 + counts1
  if ( is.matrix(reference) ) y = rowSums(1+reference)
  else y = 1 + reference
  M = log10(x/y)
  A = log10(x*y)/2
  lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
  loM = predict(lo, A)
  loM = loM - sum(loM*10^(A))/sum(10^(A))

  counts1 = counts1*10^(-loM)
  return(counts1)
}

#helper functions that extracts information from annotation objects.
annotationToX = function(annotation, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(chrLengths(genome)), 'outside')
  chr =gsub('chr','',gsub(';.*', '', as.character(annotation$Chr)))
  start = as.numeric(gsub(';.*', '', annotation$Start))
  end =  as.numeric(gsub('.*;', '', annotation$End))
  x = (end + start)/2 + prevChrL[chr]

  return(x)
}
annotationToX1X2 = function(annotation, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(chrLengths(genome)), 'outside')
  chr =gsub('chr','',gsub(';.*', '', as.character(annotation$Chr)))
  start = as.numeric(gsub(';.*', '', annotation$Start))
  end =  as.numeric(gsub('.*;', '', annotation$End))
  x1 = start + prevChrL[chr]
  x2 = end + prevChrL[chr]

  return(data.frame(x1, x2))
}
annotationToChr = function( annotation ) {
  return(gsub('chr','',gsub(';.*', '',annotation$Chr)))
}


#' a plotting function for MA plots.
#'
#' @details This plotting function does normal MA plots, with nice density visuals. Also have a range of options for things such as library normalisation and loess fits.
#'
#' @param x The x coordinates. Large x will push points towards up right.
#' @param y The y coordinates. Large y will push points towards up right.
#' @param col The base colour of the dots.
#' @param libNorm If library normalisation should be perform on x and y before plotting.
#' @param span The span of the loess fit. Higher values gives smoother fits.
#' @param medianSigma Boolean. Prints the median deviation from 0 compared to poisson standard deviation. This is a good measure of how much systematic differences there are between x and y.
#' @param loess Boolean. Plots a loess fit of the data.
#' @param verbose Boolean. Outputs some messages while plotting.
#' @param ... Remaining arguments are passed on to plot(...).
#'
#' @export
#' @examples
#' x = rpois(3000, 1:3000)
#' y = rpois(3000, 1:3000)
#' plotMA(x, y)
#' plotMA(x, y, loess=T)
#'
#' x = rpois(3000, round((1:3000)*1.2))
#' y = rpois(3000, 1:3000)
#' plotMA(x, y, loess=T)
#' 
plotMA = function(x, y, col=mcri('darkblue'), libNorm = F, span=0.2, medianSigma=T,
  xlab='A = log2(x*y)/2', ylab='M = log2(x/y)', loess=F, cex=0.6, pch=16, verbose=T, ...) {
  Nx = sum(x)
  Ny = sum(y)
  if ( libNorm ) {
    x = x/Nx
    y = y/Ny
  }
  xmin = ymin = 1
  if ( libNorm ) {
    xmin = smear/Nx
    ymin = smear/Ny
  }
  x = xmin*0.2 + noneg(xmin*0.8 + x + pmax(-0.45*xmin, pmin(0.45*xmin, rnorm(length(x), 0, xmin*0.15))))
  y = ymin*0.2 + noneg(ymin*0.8 + y + pmax(-0.45*xmin, pmin(0.45*xmin, rnorm(length(y), 0, ymin*0.15))))
  A = log2(x*y)/2
  M = log2(x/y)
  plot(A, M, cex=cex, pch=pch, xlab=xlab, ylab=ylab, yaxt='n', col=col, ...)
  segments(rep(-30, 5), c(0,-1,1, log2(c(10, 0.1))), rep(30,5), c(0,-1,1, log2(c(10, 0.1))), col=rgb(.5, .5, .5, .2), lwd=2)
    points(A, M, cex=2/3*cex, pch=pch,
           col=mcri('blue', 0.4))
    points(A, M, cex=1/2*cex, pch=pch,
           col=mcri('azure', 0.1))
    points(A, M, cex=1/3*cex, pch=pch,
           col=mcri('green', 0.02))
  axis(2, at=c(log2(0.1), -1, 0,1,log2(10)), labels=c('log2(0.1)', '-1', '0', '1', 'log2(10)'), cex.axis=1)
    
  if ( loess ) {
    if ( verbose ) catLog('Calculating loess fit...')
    lo = loess(M~A, data.frame(M, A), control=loess.control(trace.hat = 'approximate'), span=span)
    if ( verbose ) catLog('done.\n')
    As = min(A) + (0:100)/100*(max(A) - min(A))
    lines(As, predict(lo, As), col=mcri('orange'), lwd=6)
  }

  if ( medianSigma ) {
    poisErr = log2(1+2^(-A/2))
    sigma = abs(M)/poisErr
    legend('topright', paste0('<s>=', round(mean(sigma), 2)))
  }
}



#' A better version of plot
#'
#' @details This plotting function does the same as plot() in essence, but the default looks nicer than the black rings, and overplotting smoothly transitions into heatmap representation. The price is that pdfs are a bit larger than from a plot() call.
#'
#' @param x Numeric. The x coordinates.
#' @param y Numeric. The y coordinates.
#' @param add boolean. If adding the plot onto whatever is alrady there. Essentially turns plot() it into points().
#' @param showDensity Boolean. Adds a gradient through cyan and green in dense regions. Default TRUE.
#' @param ... Remaining parameters are passed to plot(...), or to points(...) if add.
#'
#' @export
#' @examples
#' x = rnorm(10000, ((1:10000)/10000)^2, 0.01)
#' y = rnorm(10000, ((1:10000)/10000)^3, 0.01)
#' plotColourScatter(x, y)
#'
plotColourScatter = function(x, y, xlab='', ylab='', col='defaultBlue', main='cor',
  add=F, cex=1, showDensity=T,...) {
  if ( main == 'cor' ) main = paste('Correlation is', signif(cor(x,y), 2))
  if ( col == 'defaultBlue' ) plotCol=mcri('darkblue')
  else if ( col == 'defaultRed' ) plotCol=mcri('darkred')
  else plotCol = col
  if ( !add ) plot(x, y, cex=cex*0.6, pch=16, xlab=xlab, ylab=ylab,
                   col=plotCol, main=main, ...)
  else points(x, y, cex=cex*0.6, pch=16, col=plotCol, ...)
  if ( showDensity & col == 'defaultBlue' ) {
    points(x, y, cex=cex*0.4, pch=16, col=mcri('blue', 0.4))
    points(x, y, cex=cex*0.3, pch=16, col=mcri('azure', 0.1))
    points(x, y, cex=cex*0.2, pch=16, col=mcri('green', 0.02))
  }
  if ( showDensity & col == 'defaultRed' ) {
    points(x, y, cex=cex*0.4, pch=16, col=mcri('red', 0.4))
    points(x, y, cex=cex*0.3, pch=16, col=mcri('orange', 0.1))
    points(x, y, cex=cex*0.2, pch=16, col=mcri('white', 0.02))
  }
}

#helper function that replaces negative values with 0.
noneg = function(x) return(ifelse(x < 0, 0, x))

#helper wrappers of colSums etc that handle non-matrices.
colsums = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (colSums(mx))
}
rowsums = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (rowSums(mx))
}
colmeans = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (colMeans(mx))
}
rowmeans = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  else (rowMeans(mx))
}

#normalise all columns to the average column size
libNorm = function(mx) {
  if ( !is.matrix(mx) ) return(mx)
  sizes = colSums(mx)
  av = mean(colSums(mx))
  ret = t(t(mx)*av/sizes)
  return(ret)
}
