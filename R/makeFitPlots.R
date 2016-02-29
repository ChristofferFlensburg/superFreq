

#' Plots results from differential coverage analysis
#'
#' @importFrom WriteXLS WriteXLS
makeFitPlots = function(fit, plotDirectory, genome, forceRedoVolcanoes=F, forceRedoDifferentRegions=F) {
  dirname = paste0(plotDirectory, '/volcanoes/')
  if ( !file.exists(dirname) ) dir.create(dirname)

  catLog('Plotting volcanoes to ', dirname, '..', sep='')
  for (col in colnames(fit$fit) ) {
    volcanoFile = paste0(dirname, col, '.jpg')
    if ( !file.exists(volcanoFile) | forceRedoVolcanoes ) {
      catLog(col, '..', sep='')
      jpeg(volcanoFile, height = 10, width = 15, res=300, units='in')
      plotVolcano(fit$fit, coef=col)
      dev.off()
    }
    volcanoFile = paste0(dirname, col, '.exon.jpg')
    if ( !file.exists(volcanoFile) | forceRedoVolcanoes ) {
      catLog(col, '..', sep='')
      jpeg(volcanoFile, height = 10, width = 15, res=300, units='in')
      try(plotVolcano(fit$exonFit, coef=col))
      dev.off()
    }
  }
  catLog('done!\n', sep='')


  differentRegionFile = paste0(plotDirectory, '/differentRegionsSamples.xls')
  if ( !file.exists(differentRegionFile) | forceRedoDifferentRegions ) {
    catLog('Writing different regions to ', differentRegionFile, '..', sep='')
    tops = list()
    for (col in colnames(fit$fit) ) {
      catLog(col, '..', sep='')
      fdr = p.adjust(fit$fit$p.value[,col], method='fdr')
      ord = order(fit$fit$XRank[,col])
      if ( length(ord) > 65534 ) {
        ord = ord[1:65534]
        catLog('outputting top 65k DE regions only, for excel...')
      }
      
      tops[[col]] = data.frame(
            'captureRegion'=rownames(fit$fit),
            'chr'=xToChr(fit$fit$x, genome=genome),
            'pos'=xToPos(fit$fit$x, genome=genome),
            'LFC'=fit$fit$coefficients[,col],
            'width'=abs(fit$fit$coefficients/fit$fit$t)[,col],
            'p-value'=fit$fit$p.value[,col],
            'FDR'=fdr,
            'correctedLFC'=fit$fit$best.guess[,col],
            'XRank'=fit$fit$XRank[,col]
            )[ord,]
    }
    names(tops) = substring(names(tops), 1, 30)
    WriteXLS('tops', differentRegionFile)
    catLog('done!\n', sep='')
  }

  differentRegionFile = paste0(plotDirectory, '/differentRegionsSamples.exons.xls')
  if ( !file.exists(differentRegionFile) | forceRedoDifferentRegions ) {
    catLog('Writing different regions to ', differentRegionFile, '..', sep='')
    tops = list()
    for (col in colnames(fit$exonFit) ) {
      catLog(col, '..', sep='')
      fdr = p.adjust(fit$exonFit$p.value[,col], method='fdr')
      ord = order(fit$exonFit$XRank[,col])
      if ( length(ord) > 65534 ) {
        ord = ord[1:65534]
        catLog('outputting top 65k DE regions only, for excel...')
      }
      
      tops[[col]] = data.frame(
            'captureRegion'=rownames(fit$exonFit),
            'chr'=xToChr(fit$exonFit$x, genome=genome),
            'pos'=xToPos(fit$exonFit$x, genome=genome),
            'LFC'=fit$exonFit$coefficients[,col],
            'width'=abs(fit$exonFit$coefficients/fit$exonFit$t)[,col],
            'p-value'=fit$exonFit$p.value[,col],
            'FDR'=fdr,
            'correctedLFC'=fit$exonFit$best.guess[,col],
            'XRank'=fit$exonFit$XRank[,col]
            )[ord,]
    }
    names(tops) = substring(names(tops), 1, 30)
    WriteXLS('tops', differentRegionFile)
    catLog('done!\n', sep='')
  }
}

#helepr functions converting from genomic coordinates to chr+bp
xToChr = function(x, genome='hg19') {
  chrL = chrLengths(genome)
  ret = rep('', length(x))
  for ( chr in names(chrL) ) {
    ret[x > 0 & x < chrL[chr]] = chr
    x = x - chrL[chr]
  }
  return(ret)
}
xToPos = function(x, genome='hg19') {
  chr = xToChr(x, genome)
  pos = x - cumsum(chrLengths(genome))[chr] + chrLengths(genome)[chr]
  return(pos)
}
