
#plots the CNV calls and clonalities, and the underlying coverage and frequencies that the calls are based on.
makeCNVplots = function(cnvs, plotDirectory, genome='hg19', plotPDF=F, forceRedoCNVplots=F) {
  CNVplotDirectory = paste0(plotDirectory, '/CNV/')
  if ( !file.exists(CNVplotDirectory) ) dir.create(CNVplotDirectory)
  for ( name in names(cnvs) ) {
    dirname = paste0(CNVplotDirectory, name)
    if ( !file.exists(dirname) ) dir.create(dirname)

    #plot individual capture regions
    if ( plotPDF )
      filename = paste0(dirname, '/captureRegions.pdf')
    else
      filename = paste0(dirname, '/captureRegions.png')
    if ( !file.exists(filename) | forceRedoCNVplots ) {
      catLog('Plotting CNVs to ', dirname, '.\n', sep='')
      if ( plotPDF )
        pdf(filename, width=20, height=10)
      else
        png(filename, width=20, height=10, res=300, units='in')
      plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome)
      dev.off()
    }

    #plot merged regions
    if ( plotPDF )
      filename = paste0(dirname, '/CNVcalls.pdf')
    else
      filename = paste0(dirname, '/CNVcalls.png')
    if ( !file.exists(filename) | forceRedoCNVplots ) {
      if ( plotPDF )
        pdf(filename, width=20, height=10)
      else
        png(filename, width=20, height=10, res=300, units='in')
      plotCR(cnvs[[name]]$clusters, genome=genome)
      dev.off()
    }

    

    #plot both
    if ( plotPDF )
      filename = paste0(dirname, '/combined.pdf')
    else
      filename = paste0(dirname, '/combined.png')
    if ( !file.exists(filename) | forceRedoCNVplots ) {
      if ( plotPDF )
        pdf(filename, width=20, height=10)
      else
        png(filename, width=20, height=10, res=300, units='in')
      plotCR(cnvs[[name]]$clusters, errorBars=F, genome=genome, alpha=0)
      plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome, alpha=0.1, add=T, moveHet=F)
      plotCR(cnvs[[name]]$clusters, errorBars=T, genome=genome, add=T)
      dev.off()
    }
    if ( plotPDF )
      outerFilename = paste0(CNVplotDirectory, '/', name, '.pdf')
    else
      outerFilename = paste0(CNVplotDirectory, '/', name, '.png')
    if ( !file.exists(outerFilename) | forceRedoCNVplots ) {
      system(paste0('cp ', gsub(' ', '\\\\ ',filename), ' ', gsub(' ', '\\\\ ',outerFilename)), intern=T)
    }

    for( chr in names(chrLengths(genome)) ) {
      if ( plotPDF )
        filename = paste0(dirname, '/chr', chr, '.pdf')
      else
        filename = paste0(dirname, '/chr', chr, '.png')
      if ( !file.exists(filename) | forceRedoCNVplots ) {
        if ( plotPDF )
          pdf(filename, width=20, height=10)
        else
          png(filename, width=20, height=10, res=300, units='in')
        plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome, chr=chr, alpha=0)
        plotCR(cnvs[[name]]$CR, errorBars=F, genome=genome, chr=chr, alpha=0.3, add=T)
        plotCR(cnvs[[name]]$clusters, errorBars=T, genome=genome, chr=chr, add=T)
        dev.off()
      }
    }
  }
}





#' plots a copy number profile
#'
#' @param cR data.frame. The "clusters" or "CR" from the data.
#' @param showClonality logical. if the clonality panel is shown, if data available. Default TRUE.
#' @param errorBars logical. If errorbars are plotted. Default TRUE.
#' @param chr character. Which chromosome to plot. Default 'all', which plots the entire genome.
#' @param genome character. The genome assembly. Default 'hg19'.
#' @param alpha numeric. The opacity of the point. Default 1.
#' @param add logical. If the data should be plotted on top of whatever is already there. Default FALSE.
#' @param moveHet logical. If the SNP data should snap to f=0.5 if not significantly different. Default TRUE.
#' @param pt.cex numerical. Scaling factor for the size of the points. Default 1.
#' @param setMargins logical. If the margins should be removed to use the entire plottable area. Default TRUE.
#' @param fullFrequency logical. If the SNPs should be copied back up to 1-f as well. Generally not a good idea, and potentially misleading, but can be useful for people that are not used to see frequencies mirrored down to 0-0.5.
#' @param colourDeviation logical. If the points should be coloured based on deviation from diploid. Defaulty TRUE.
#' @param forceCol colour. The colour all the points will be plotted in. Can be useful if overplotting several samples with add=T. Overrides colourDeviation if not NA. Default NA.
#' @param plotCall logical. If the copy number call should be added to the bottom of the SNP panel. Default TRUE.
#' @param ...  remaining arguments are passed to base plot(...) if add=F, otherwise ignored. 
#'
#' @details This function plots copy number information from superFreq, either on a by-gene level or by-segment.
#'
#'
#' @export
#'
plotCR = function(cR, showClonality=T, errorBars=T, chr='all', genome='hg19', alpha=1, add=F, moveHet=T, pt.cex=1, setMargins=T, fullFrequency=F, colourDeviation=T, forceCol=NA, plotCall=T, plotArrows=F, smallPlot=F, lwd=1, sep.lwd=5, sideSpace=NULL, noSex=F, ...) {
  showClonality = showClonality & 'subclonality' %in% names(cR)
  if ( nrow(cR) == 0 ) return()
  if ( noSex ) cR = cR[!(xToChr(cR$x1, genome=genome) %in% c('X', 'Y')),]
  if ( chr != 'all' ) {
    chrL = chrLengths(genome=genome)
    xlim = c(cumsum(chrL)[chr]-chrL[chr], cumsum(chrL)[chr])
    cR = cR[cR$x1 >= xlim[1] & cR$x2 <= xlim[2],]
    if ( nrow(cR) == 0 ) return()
  }
  else xlim=c(min(cR$x1), max(cR$x2))
  xlimMar = xlim
  if ( is.null(sideSpace) ) {
    sideSpace = 0.05
    if ( smallPlot ) sideSpace = 0.11
  }
  xlimMar[1] = xlim[1] - (xlim[2]-xlim[1])*sideSpace
  xlimMar[2] = xlim[2] + (xlim[2]-xlim[1])*sideSpace

  ymax = 1.1
  if ( smallPlot ) ymax = 1.2
  ylim = c(-2.3, ymax)
  if ( !showClonality ) ylim = c(-1.1, ymax)
  if ( !add ) {
    if ( setMargins ) {
      par(oma=rep(0, 4))
      par(mar=c(0, 0, 0, 0))
    }
    plot(0, xlim=xlimMar, ylim=ylim, type='n', xlab='', ylab='', axes=F, ...)
    if ( fullFrequency ) {
      lineA = 0.3
      labelA = 1
      lineY = c(0.1,0.6,0.6+log2(1.5)/2,1.1, -0.1,-1.1+2/3, -0.6,-1.1+1/3,-1.1)
      lineCol = c(rgb(0,0,0,lineA), rgb(0,0.8,0,lineA), rgb(1,0.5,0,lineA), rgb(1,0,0,lineA),
        rgb(0,0,0,lineA), rgb(1,0.5,0,lineA), rgb(1,0,0,lineA), rgb(1,0.5,0,lineA), rgb(0,0,0, lineA))
      labelCol = c(rgb(0,0,0,labelA), rgb(0,0.8,0,labelA), rgb(1,0.5,0,labelA), rgb(1,0,0,labelA),
        rgb(0,0,0,labelA), rgb(1,0.5,0,labelA), rgb(1,0,0,labelA), rgb(1,0.5,0,labelA), rgb(0,0,0, labelA))
      lineName = c('1 chr', '2 chr', '3 chr', '4 chr', '100%', '67%', '50%', '33%', '0%')
      lowerLabel = 'allele frequency'
      lowerTicks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    }
    else {
      lineA = 0.3
      labelA = 1
      lineY = c(0.1,0.6,0.6+log2(1.5)/2,1.1, -0.1,-1.1+2/3, -0.6, -1.1)
      lineCol = c(rgb(0,0,0,lineA), rgb(0,0.8,0,lineA), rgb(1,0.5,0,lineA), rgb(1,0,0,lineA),
        rgb(0,0.8,0,lineA), rgb(1,0.5,0,lineA), rgb(1,0,0,lineA), rgb(0,0,0, lineA))
      labelCol = c(rgb(0,0,0,labelA), rgb(0,0.8,0,labelA), rgb(1,0.5,0,labelA), rgb(1,0,0,labelA),
        rgb(0,0.8,0,labelA), rgb(1,0.5,0,labelA), rgb(1,0,0,labelA), rgb(0,0,0, labelA))
      lineName = c('1 chr', '2 chr', '3 chr', '4 chr', '50%', '33%', '25%', '0%')
      lowerLabel = 'SNP MAF'
      lowerTicks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
    }
    segments(xlim[1], lineY, xlim[2], lineY, lwd=3, col=lineCol)
    rightShift = 0.005
    if ( smallPlot ) rightShift = 0.02
    text(xlim[2] + (xlim[2]-xlim[1])*rightShift, lineY, lineName, adj=0, cex=1, col=labelCol)
    if ( smallPlot )
      text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.02, 0.6, srt=90, 'LFC vs normals', cex=0.8)
    else
      text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.02, 0.6, srt=90, 'coverage LFC vs normals', cex=0.8)
    text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.02, -0.6, srt=90, lowerLabel, cex=0.8)
    text(xlimMar[1] - (xlimMar[2]-xlimMar[1])*0.02, -1.8, srt=90, 'clonality', cex=0.8)
    if ( smallPlot ) addChromosomeLines(ylim=c(-2.3, 1.18), col=mcri('green', 0.6), lwd=1, genome=genome, onlyNumbers=noSex)
    else addChromosomeLines(ylim=c(-2.3, 1.15), col=mcri('green', 0.6), lwd=1, genome=genome, onlyNumbers=noSex)
    segments(2*xlim[1]-xlim[2], c(0, -1.2), 2*xlim[2]-xlim[1], c(0, -1.2), lwd=sep.lwd)
    axis(side=2, at=0.1 + (0:4)/4, labels=c(-1, -0.5, 0, 0.5, 1), pos=xlim[1], padj=0.75)
    axis(side=2, at=-1.1 + (0:5)/5, labels=lowerTicks, pos=xlim[1], padj=0.75)
    axis(side=2, at=-2.3 + (0:5)/5, labels=c(0, 0.2, 0.4, 0.6, 0.8, 1), pos=xlim[1], padj=0.75)
  }
  
  x = (cR$x1 + cR$x2)/2
  y = 0.6+cR$M/2
  col = rgb(pmin(1,pmax(0,2*cR$M))*colourDeviation,0,pmin(1,pmax(0,-2*cR$M))*colourDeviation, alpha)
  if ( !is.na(forceCol)[1] ) col = forceCol 
  width = sqrt(cR$width^2+getSystematicVariance()^2)
  lfcCex = pt.cex*pmin(3, sqrt(0.1/width))
  points(x, pmax(0, y), cex=lfcCex, col=col, pch=16)
  if ( errorBars ) segments(x, pmax(0,y-width/2), x, pmax(0, y+width/2), lwd=lfcCex*lwd, col=col)
  segments(cR$x1, pmax(0,y), cR$x2, pmax(0, y), lwd=lfcCex*lwd, col=col)
  tooHigh = y > 1.2
  tooLow = y < 0
  if ( any(tooHigh) & plotArrows ) {
    highCol = rgb(1*colourDeviation,0,0,alpha)
    if ( !is.na(forceCol)[1] ) highCol = forceCol 
    arrows(x[tooHigh], 1.12, x[tooHigh], 1.18, lwd=2, length=0.1, col=highCol)
    text(x[tooHigh], 1.09, round(cR$M[tooHigh],2), cex=0.8, col=highCol)
  }
  if ( any(tooLow) & plotArrows ) {
    lowCol = rgb(0,0,1*colourDeviation,alpha)
    if ( !is.na(forceCol)[1] ) lowCol = forceCol 
    arrows(x[tooLow], 0.11, x[tooLow], 0.05, lwd=2, length=0.1, col=lowCol)
    text(x[tooLow], 0.14, round(cR$M[tooLow],2), cex=0.8, col=lowCol)
  }

  xf = (cR$x1+cR$x2)[cR$cov > 0]/2
  cRf = cR[cR$cov > 0,]
  if ( nrow(cRf) > 0 ) {
    f = refUnbias(cRf$var/cRf$cov)
    if ( 'f' %in% colnames(cR) ) f = cRf$f
    yf = -1.1 + 2*f
    #to protect from 0 variant reads giving 0 error below, we pmax with a Poisson, for getting 0 from
    #a non-zero expected alt count.
    #getting 0 from Poisson with lambda=2 is p-value 86%, which roughly corresponds to the CI of our error.
    if ( 'ferr' %in% names(cR) ) ferr = cRf$ferr
    else ferr = 2*sqrt(pmax(2,cRf$var))/cRf$cov
    fCex = pt.cex*pmax(0.2, pmin(3, sqrt(0.025/ferr)))
    #plot the average MAFs, opaqueness from how likely a non-het is.
    pHet = cRf$pHet
    if ( 'postHet' %in% colnames(cRf) ) pHet = cRf$postHet
    if ( 'stat' %in% colnames(cRf) ) pHet[cRf$stat > 0] = 1
    if ( 'call' %in% colnames(cRf) ) {
      pHet = rep(0, length(pHet))
      pHet[cRf$call %in% c('AB', 'CL', 'AABB')] = 1
    }
    if ( moveHet ) {
      col = rgb(pmax(0,pmin(1, 4*(0.5-f)))*colourDeviation,0,0, (1-pHet)*alpha)
      if ( !is.na(forceCol)[1] ) col = ifelse(pHet*alpha <= 0.5, forceCol, rgb(0,0,0,0))
    }
    else {
      col = rgb(pmax(0,pmin(1, 4*(0.5-f)))*colourDeviation,0,0, alpha)
      if ( !is.na(forceCol)[1] ) col = forceCol 
    }
    if ( fullFrequency ) {
      xf = rep(xf, 2)
      yf = c(-1.1 + f, -0.1 - f)
      col = rep(col, 2)
      ferr = rep(ferr, 2)
    }
    points(xf, yf, cex=fCex, pch=16, col=col)
    if ( errorBars ) segments(xf, pmin(-0.1, yf+ferr), xf, pmax(-1.1, yf-ferr), col=col, lwd=fCex*lwd)
    segments(cRf$x1, yf, cRf$x2, yf, lwd=fCex*lwd, col=col)
    #plot f=0.5, opaqueness from how likely a het is.
    if ( moveHet ) {
      col = rgb(0,0,0, pHet*alpha)
      if ( !is.na(forceCol)[1] ) col = ifelse(pHet*alpha > 0.5, forceCol, rgb(0,0,0,0))
      yf = rep(-0.1, length(xf))
      if ( fullFrequency ) {
        xf = (cR$x1+cR$x2)[cR$cov > 0]/2
        yf = rep(-0.6, length(xf))
        ferr = ferr[1:length(yf)]
      }
      points(xf, yf, cex=fCex, pch=16, col=col)
      if ( errorBars ) segments(xf, pmin(-0.1, yf+ferr), xf, pmax(-1.1, yf-ferr), col=col, lwd=fCex*lwd)
      segments(cRf$x1, yf, cRf$x2, yf, lwd=fCex*lwd, col=col)
    }
  }

  if ( 'call' %in% colnames(cR) & plotCall ) {
    called = which(cR$call != 'AB')
    if ( length(called) > 0 ) {
      y = rep(c(-1.05, -0.95, -0.85), length(called))[1:length(called)]
      calls = gsub(' ', '', shortenCalls(paste0(' ', cR$call[called])))
      text(x[called], y-0.05, calls, cex=0.9, col=mcri('green', alpha))
    }
  }

  if ( showClonality ) {
    subcloneMax = sort(unique(cR$subclonality + cR$subclonalityError), decreasing=T)
    subclone = sort(unique(cR$subclonality), decreasing=T)
    subcloneMin = sort(unique(cR$subclonality - cR$subclonalityError), decreasing=T)
    subcloneSigma = sapply(subclone, function(subF) mean(((cR$clonality-subF)/cR$clonalityError)[cR$subclonality == subF]))
    if ( length(subclone) == 1 ) subcloneCol = subcloneColBG = rgb(0,0,0,alpha)
    else {
      subcloneCol = c(rgb(0,0,0,alpha), randomCol(1:(length(subclone)-1), a=alpha, noBlack=T))
      subcloneColBG = c(rgb(0,0,0,alpha), randomCols(1:(length(subclone)-1),
        pmax(0.02, pmin(0.25, (0.01/pmax(subcloneSigma,1)/(subcloneMax-subcloneMin))[-1]))*alpha, noBlack=T))
    }
    col = sapply(cR$subclonality, function(sub) subcloneCol[which(subclone == sub)])
    rect(xlim[1], subcloneMin - 2.3, xlim[2], pmin(-1.3, subcloneMax -2.3),
         col=subcloneColBG, border=F)
    y = -2.3 + cR$clonality
    error = cR$clonalityError
    clonCex = pt.cex*pmax(0.2, pmin(2, sqrt(0.05/error)))
    points(x, y, pch=16, cex=clonCex, col=col)
    segments(x, y-error, x, pmin(-1.3, y+error), lwd=clonCex*lwd, col=col)
    segments(cR$x1, y, cR$x2, y, lwd=clonCex*lwd, col=col)
  }
}

#helper function that adds vertical chromsome lines to a plot on the genomic x coordinate.
addChromosomeLines = function(ylim = c(0, 1), col = 'red', cex=1, genome='hg19', onlyNumbers=F,
							  stagger=0, staggerFrom=1, excludeChromosomes=c(), ...) {
  chrL = chrLengths(genome=genome)
  if ( onlyNumbers ) chrL = chrL[grepl('[0-9]', names(chrL))]
  chrL = chrL[!(names(chrL) %in% excludeChromosomes)]
  lim = c(0, cumsum(chrL))
  segments(lim, ylim[1], lim, ylim[2], col=col, cex=cex, ...)
  ys = rep(ylim[2], length(chrL))
  if ( stagger != 0 ) {
    is = staggerFrom:length(chrL)
    ys[is] = ys[is] + rep(c(stagger,0),length(is))[1:length(is)]
  }
  text((lim[-1] + lim[-length(lim)])/2, ys, names(chrL), col=col, cex=cex)
}

#Helper function that returns a colour (or line type) depending on the provided index. Cyclic.
#colours (with the no-black option) and line types have cycles of 6 and 7, so shouldn't repeat until index 43.
randomLtys = function( i ) {
  ltys = c(1,3,5,4,2,6)
  return(ltys[1+(i-1)%%length(ltys)])
}
randomCol = function( i, a=1, noBlack=F ) {
  cols = c(rgb(1,0,0,a), rgb(0,0,1,a), rgb(0, 0.7,0,a), rgb(0,0,0,a), rgb(0.9, 0.6,0,a), rgb(0.9,0,0.9,a), rgb(0,0.8,0.7,a),
    rgb(0.6,0.6,0.6,a))
  if ( noBlack ) cols = cols[-c(4,8)]
  return(cols[1+(i-1)%%length(cols)])
}
randomCols = function(i, a=1, noBlack=F) {
  if ( length(a)==1 & length(i) > 1 ) a = rep(a, length(i))
  if ( length(i) != length(a) ) cat('randomCols wants colour index and alpha of same length! >:(\n')
  return(sapply(1:length(i), function(I) randomCol(i[I], a[I], noBlack=noBlack)))
}


#helper function that restores default margins
resetMargins = function() {
  par(oma=c(0,0,0,0))
  par(mar=c(4,4,3,1))
}


#makes the multi-sample CNA summary heatmap across the genome.
#WIP at this point.
plotCNAheatmap = function(clusters, genome, onlyNumbers=F, add=F, yShift=0, legend=T, xlim=NULL) {
  superFreq:::resetMargins()
  chrL = chrLengths(genome)
  if ( is.null(xlim) ) xlim=c(0, 1)*sum(chrL)
  xlimPlot = c(xlim[1] - 0.1*(xlim[2]-xlim[1]), xlim[2])
  if ( !add ) plot(1, type='n', xlim=xlimPlot, ylim=c(0.5, 1+length(clusters)), xlab='', ylab='', xaxt='n', yaxt='n', frame.plot=F)
  samples = names(clusters)
  for ( sample in samples ) {
    i = which(sample == samples)
    cluster = clusters[[sample]]$clusters
    cluster = cluster[cluster$x2 > xlim[1] & cluster$x1 < xlim[2],]
    if ( onlyNumbers ) cluster = cluster[!(xToChr(cluster$x1, genome=genome) %in% c('X', 'Y')),]
    isCNN = cluster$call %in% c('AA', 'AA?', 'AA??')
    isAB = cluster$call == 'AB'
    red=superFreq:::noneg(pmin(cluster$M,1))*(!isCNN)*(!isAB)
    darkred = cluster$M > 1 | nchar(gsub('?', '', cluster$call)) > 4
    blue=superFreq:::noneg(pmin(-cluster$M,1))*(!isCNN)*(!isAB)
    green = isCNN*(0.5-cluster$f)*2*(!isAB)
    green = ifelse(is.na(green), 0, green)
    col = ifelse(darkred, mcri('darkred'), ifelse(red > 0, mcri('red', red), ifelse(blue > 0, mcri('blue', blue), mcri('green', green))))
    width = ifelse(red+green+blue > 0.3, 1, 0)
    #dont plot outside xlim
    x1 = pmax(cluster$x1, xlim[1])
    x2 = pmin(cluster$x2, xlim[2])

    y = yShift + i
    rect(x1, y-width/2, x2, y+width/2, col=col, border=F)
    text(xlim[1] -0.004*(xlim[2]-xlim[1]), y, sample, adj=1, cex=0.7)
  }
  if ( !add ) superFreq:::addChromosomeLines(ylim=c(0.5, length(clusters)+1), col=mcri('grey'), genome=genome, onlyNumbers=onlyNumbers)
  if ( !add & legend ) legend('bottomright', c('amp', 'gain', 'loss', 'LOH'), pch=15, pt.cex=1.7, col=mcri(c('darkred', 'red', 'blue', 'green')), bg='white', cex=0.8)
}

#makes summary plots of CNAs across all individuals in the batch.
makeCNAheatmap = function(clusters, plotDirectory, genome) {
  summaryDir = paste0(plotDirectory, '/summary')
  superFreq:::ensureDirectoryExists(summaryDir)
  plotFile = paste0(summaryDir, '/CNAheatmap.pdf')
  if ( file.exists(plotFile) ) return()
  pdf(plotFile, width=15, height=10)
  plotCNAheatmap(clusters=clusters, genome=genome)
  dev.off()
}


#merges and makes CNA heatmaps across all the individuals of a batch.
#requires the batch to be run with split=T (Which you should anyway).
plotCNAbatchHeatmap = function(Rdirectory, metaDataFile, genome, onlyNumbers=F, xlim=NULL, excludeIndividuals=c(), excludeSamples=c(), cpus=1) {
  clusterList = superFreq:::loadClusterList(Rdirectory=Rdirectory, metaDataFile=metaDataFile, excludeIndividuals=excludeIndividuals,
                                excludeSamples=excludeSamples, cpus=cpus)

  chrL = chrLengths(genome)
  ymax = 1+sum(sapply(clusterList, length)) + 0.5*(length(clusterList)-1)
  if ( is.null(xlim) ) xlim=c(0, 1)*sum(chrL)
  xlimPlot = c(xlim[1] - 0.1*(xlim[2]-xlim[1]), xlim[2])
  plot(1, type='n', xlim=xlimPlot, ylim=c(0.5, ymax), xlab='', ylab='', xaxt='n', yaxt='n', frame.plot=F)
  for ( i in 1:length(clusterList) ) {
    yShift=c(0,cumsum(sapply(clusterList, length)))[i] + 0.5*(i-1)
    superFreq:::plotCNAheatmap(clusterList[[i]], genome=genome, onlyNumbers=onlyNumbers, add=T, yShift=yShift, xlim=xlim)
    #if many samples, thinner separating lines or no lines at all.
    if ( i > 1 & length(clusterList) < 100 )
      segments(xlimPlot[1], yShift + 0.25, xlimPlot[2], yShift + 0.25, lwd=min(1, 30/length(clusterList)), col='grey')
  }
  superFreq:::addChromosomeLines(ylim=c(0.5, ymax*1.01), col=mcri('grey'), genome=genome, onlyNumbers=onlyNumbers)
  legend('bottomright', c('amp', 'gain', 'loss', 'LOH'), pch=15, pt.cex=1.7, col=mcri(c('darkred', 'red', 'blue', 'green')), bg='white', cex=0.8)
}

#makes summary plots of CNAs across all individuals in the batch.
makeCNAbatchHeatmap = function(Rdirectory, plotDirectory, metaDataFile, genome, excludeIndividuals=c(), excludeSamples=c(), cpus=1) {
  cohortDir = paste0(plotDirectory, '/cohortWide')
  superFreq:::ensureDirectoryExists(cohortDir)
  plotFile = paste0(cohortDir, '/CNAsummary.pdf')
  md = importSampleMetaData(metaDataFile)
  md = md[!(md$INDIVIDUAL %in% excludeIndividuals) & !(md$NAME %in% excludeSamples),]
  pdf(plotFile, width=15, height=max(7, nrow(md)/5))
  superFreq:::plotCNAbatchHeatmap(Rdirectory=Rdirectory, metaDataFile=metaDataFile, genome=genome,
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  dev.off()
}
