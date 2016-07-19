
#takes the subclone evolution from the stories and plots rivers, showing which events are aprt of which subclone.
makeRiverPlots = function(stories, variants, genome='hg19', cpus=1, plotDirectory, forceRedo=F) {
  if ( length(stories) == 0 ) return
  riverDirectory = paste0(plotDirectory, '/rivers/')
  if ( !file.exists(riverDirectory) ) dir.create(riverDirectory)
  for ( ts in names(stories) ) {
    file = paste0(riverDirectory, ts, '-river.pdf')
    if ( file.exists(file) & !forceRedo ) next
    catLog('Making riverplot for ', ts, '..', sep='')
    patVar = variants
    patVar$variants = patVar$variants[colnames(stories[[ts]]$consistentClusters$cloneStories$stories)]

    pdf(file, width=15, height=10)
    cloneCols =
      plotRiver(cloneTree=stories[[ts]]$consistentClusters$cloneTree, cloneStories=stories[[ts]]$consistentClusters$cloneStories,
                storyList=stories[[ts]]$consistentClusters$storyList, allStories=stories[[ts]]$allConsistent,
                variants=patVar, genome=genome)
    plotStories(stories[[ts]]$consistentClusters$cloneStories, patVar, col=cloneCols, genome=genome)
    heatmapStories(stories[[ts]]$allConsistent, stories[[ts]]$consistentClusters$storyList,
                   patVar, col=cloneCols, genome=genome)
    for ( subclone in names(stories[[ts]]$consistentClusters$storyList) ) {
      i = which(names(stories[[ts]]$consistentClusters$storyList) == subclone)
      plotStories(stories[[ts]]$allConsistent[stories[[ts]]$consistentClusters$storyList[[subclone]],],
                  patVar, alpha=0.2, genome=genome)
      plotStories(stories[[ts]]$consistentClusters$cloneStories[i,], patVar, add=T,
                  col=cloneCols[i], genome=genome)
    }
    dev.off()

    if ( length(stories[[ts]]$clusters$storyList) > length(stories[[ts]]$consistentClusters$storyList) ) {
      removed = !(names(stories[[ts]]$clusters$storyList) %in% names(stories[[ts]]$consistentClusters$storyList))
      catLog('showing ', sum(removed), ' dodgy clones..', sep='')
      removedFile = paste0(riverDirectory, ts, '-river-DODGY.pdf')
      pdf(removedFile, width=15, height=10)
      cloneCols = plotRiver(stories[[ts]]$cloneTree, stories[[ts]]$clusters$cloneStories,
                            stories[[ts]]$clusters$storyList, stories[[ts]]$all,
                            variants=patVar, genome=genome)
      plotStories(stories[[ts]]$clusters$cloneStories, patVar, col=cloneCols, genome=genome)
      heatmapStories(stories[[ts]]$all, stories[[ts]]$clusters$storyList,
                     patVar, col=cloneCols, genome=genome)
      for ( subclone in names(stories[[ts]]$clusters$storyList) ) {
        i = which(names(stories[[ts]]$clusters$storyList) == subclone)
        plotStories(stories[[ts]]$all[stories[[ts]]$clusters$storyList[[subclone]],],
                    patVar, alpha=0.2, genome=genome)
        plotStories(stories[[ts]]$clusters$cloneStories[i,], patVar, add=T,
                    col=cloneCols[i], genome=genome)
      }
      dev.off()
    }
    catLog('done.\n')
    
    catLog('Outputting data on stories ', ts, '..', sep='')
    excelFile = paste0(riverDirectory, ts, '-river.xls')
    output = do.call(rbind, lapply(stories[[ts]]$consistentClusters$storyList,
      function(sL) stories[[ts]]$allConsistent[sL,]))
    clone = do.call(c, lapply(1:length(stories[[ts]]$consistentClusters$storyList),
      function(i) rep(names(stories[[ts]]$consistentClusters$storyList)[i], length(stories[[ts]]$consistentClusters$storyList[[i]]))))
    output$clone = clone
    chr = xToChr(output$x1, genome)
    start = xToPos(output$x1, genome)
    end = xToPos(output$x2, genome)
    label = storyToLabel(output, patVar, genome, maxLength=100)$label
    names(start) = names(end) = names(chr) = names(label) = rownames(output)
    clonality = as.data.frame(output$stories)
    error = as.data.frame(output$errors)
    names(chr) = names(start) = names(end) = names(label) = rownames(clonality)
    output = cbind(chr = chr, start=start, end=end, name=label, clone=output$clone, clonality=clonality, error=error)
    output = addAnnotationToOutput(output, patVar, genome=genome)
    WriteXLS('output', excelFile)

    if ( length(stories[[ts]]$clusters$storyList) > length(stories[[ts]]$consistentClusters$storyList) ) {
      removed = !(names(stories[[ts]]$clusters$storyList) %in% names(stories[[ts]]$consistentClusters$storyList))
      catLog('adding dodgy clones..', sep='')
      excelFile = paste0(riverDirectory, ts, '-river-DODGY.xls')
      output = do.call(rbind, lapply(stories[[ts]]$clusters$storyList,
        function(sL) stories[[ts]]$all[sL,]))
      clone = do.call(c, lapply(1:length(stories[[ts]]$clusters$storyList),
        function(i) rep(names(stories[[ts]]$clusters$storyList)[i], length(stories[[ts]]$clusters$storyList[[i]]))))
      output$clone = clone
      chr = xToChr(output$x1, genome)
      start = xToPos(output$x1, genome)
      end = xToPos(output$x2, genome)
      label = storyToLabel(output, patVar, genome, maxLength=100)$label
      names(start) = names(end) = names(chr) = names(label) = rownames(output)
      clonality = as.data.frame(output$stories)
      error = as.data.frame(output$errors)
      output = cbind(chr = chr, start=start, end=end, name=label, clone=output$clone, clonality=clonality, error=error)
      output = addAnnotationToOutput(output, patVar, genome=genome)
      WriteXLS('output', excelFile)
    }

    catLog('done.\n')
  }
}

#helper function converting internal event names to informative labels for the plots.
storyToLabel = function(stories, variants, genome, maxLength=30) {
  call = stories$call
  isSNP = grepl('[0-9]', call)
  label = rep('', nrow(stories))
  font = rep(1, nrow(stories))
  colour = rep('black', nrow(stories))
  severity = rep(100, nrow(stories))

  q = variants$variants[[1]]
  q = q[q$x %in% stories$x1[isSNP],]  
  gene = q[stories$call,]$inGene[isSNP]
  if ( 'severity' %in% names(variants$variants[[1]]) & any(isSNP) ) {
    severityMx = sapply(variants$variants, function(q) ifelse(is.na(q[call[isSNP],]$severity), 100, q[call[isSNP],]$severity))
    if ( sum(isSNP) == 1 ) severityMx = matrix(severityMx, nrow=sum(isSNP))
    mostSevere = ceiling(apply(severityMx, 1, min))
    severity[isSNP] = mostSevere
    type = sapply(mostSevere, severityToType)
    type[type=='unknown'] = ''
    label[isSNP] = substr(paste0(gene, ' (', xToChr(stories$x1[isSNP], genome=genome), ') ', type), 1, maxLength)
  }
  else
    label[isSNP] = substr(paste0(gene, ' (', xToChr(stories$x1[isSNP], genome=genome), ') '), 1, maxLength)

  if ( 'isCosmicCensus' %in% names(variants$variants[[1]]) & any(isSNP) ) {
    cosmicMx = sapply(variants$variants, function(q) q[call[isSNP],]$isCosmicCensus)
    if ( sum(isSNP) == 1 ) cosmicMx = matrix(cosmicMx, nrow=sum(isSNP))
    isCensus = apply(cosmicMx, 1, any)
    colour[isSNP] = ifelse(isCensus, mcri('green'), 'black')
    font[isSNP] = ifelse(isCensus, 2, 1)
  }

  dist = stories$x2-stories$x1
  distText = ifelse(dist >= 1e6, paste0(round(dist/1e6), 'Mbp '), ifelse(dist >= 1e3, paste0(round(dist/1e3), 'kbp '), paste0(dist, 'bp ')))
  label[!isSNP & !is.na(dist)] = paste0(distText, call, ' (', xToChr(stories$x1, genome=genome), ')')[!isSNP & !is.na(dist)]
  label[!isSNP & !is.na(dist)] = shortenCalls(label[!isSNP & !is.na(dist)])
  font[!isSNP & !is.na(dist)] = ifelse(grepl('[0-9]AB', label[!isSNP & !is.na(dist)]), 4, 3)
  severity[!isSNP & !is.na(dist)] = 0
  #label[!isSNP & is.na(dist)] = gsub('^clone$', 'clone.0', make.names(call[!isSNP & is.na(dist)], unique=T))
  label[!isSNP & is.na(dist)] = paste0('clone.', rownames(stories)[!isSNP & is.na(dist)])
  font[!isSNP & is.na(dist)] = 4
  severity[!isSNP & is.na(dist)] = -1
  return(data.frame(label=label, colour=colour, font=font, severity=severity, stringsAsFactors=F))
}

addAnnotationToOutput = function(output, variants, genome='hg19') {
  isSNV = grepl('^[0-9]', rownames(output))
  rows = rownames(output)[isSNV]
  if ( length(rows) == 0 ) return(output)
  severityMx = sapply(variants$variants, function(q) ifelse(is.na(q[rows,]$severity), 110, q[rows,]$severity))
  if ( is.vector(severityMx) ) severityMx = matrix(severityMx, nrow=length(rows))
  mostSevere = apply(severityMx, 1, function(severities) which(severities <= 1.0001*min(severities))[1])
  columns = unique(c('severity', 'type', moreVEPnames(genome=genome)))
  for ( column in columns ) {
    if ( !(column %in% colnames(variants$variants[[1]])) ) next
    mx = sapply(variants$variants, function(q) q[rows,][,column])
    if ( is.vector(mx) ) mx = matrix(mx, nrow=length(rows))
    columnValues = mx[cbind(1:nrow(mx), mostSevere)]
    output[[column]] = rep('', nrow(output))
    output[rows,column] = columnValues
  }
  return(output)
}


#main plotting function
plotRiver = function(cloneTree, cloneStories, storyList, allStories, variants, genome='hg19', normalise=T, xlim='default', ylim='default', labels=T, setPar=T, sampleOrder='default', excludeClones=c(), markDodgy=T, colourPool = c()) {
  variants$variants = variants$variants[colnames(cloneStories$stories)]
  
  if ( length(colourPool) == 0 )
    colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))

  if ( length(sampleOrder) == 1 && sampleOrder == 'default' ) sampleOrder=colnames(cloneStories$stories)
  cloneStories$stories = cloneStories$stories[,sampleOrder,drop=F]
  cloneStories$errors = cloneStories$errors[,sampleOrder,drop=F]

  if ( length(excludeClones) > 0 ) {
    cloneTree = excludeClonesFromTree(cloneTree, excludeClones)
    cloneStories = cloneStories[!(rownames(cloneStories) %in% excludeClones), ]
    storyList = storyList[!(names(storyList) %in% excludeClones)]
  }
  
  dodgyness = getDodgyness(storyList, cloneStories)
  if ( !markDodgy ) dodgyness = 0*dodgyness
  
  cloneLabels = lapply(storyList, function(rows) storyToLabel(allStories[rows,], variants, genome=genome, maxLength=20))
  names(cloneLabels) = rownames(cloneStories)
  stories = abs(cloneStories$stories)
  rownames(stories) = rownames(cloneStories)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(abs(stories[names(cloneTree),i])))
  leadingNormal = F
  if ( any(purity==0) ) {
    leadingNormal = T
    purity[purity==0] = 1
  }
  stories[stories < cloneStories$errors*2 & stories < 0.5] = 0
  if ( normalise ) stories = t(t(stories)/purity)
  if ( !leadingNormal ) {
    stories = cbind(rep(0, nrow(stories)), stories)
    colnames(stories)[1] = 'germline'
    if ( any(rownames(stories)=='germline') ) stories['germline','germline'] = 1
  }
  x = 1:ncol(stories)

  if ( setPar ) {
    par(oma=rep(0, 4))
    par(mar=c(0, 4, 0, 0))
  }
  if ( !labels & xlim[1] == 'default' ) xlim = c(1, max(x))
  if ( labels & xlim[1] == 'default' ) xlim = c(1, max(x)+ceiling(nrow(cloneStories)/2)*1.5)
  if ( ylim[1] == 'default' ) ylim = c(-0.02,1)
  plot(1, type='n', xlim=xlim, ylim=ylim, xaxt='n', frame.plot=F,
       ylab='clonality', xlab = '')
  cloneCols = addSubclone(cloneTree, stories, ylims = matrix(rep(c(0,1), ncol(stories)), nrow=2), dodgyness, colourPool=colourPool, margin=0.02)$usedCols
  for (i in 1:nrow(cloneStories)) {
    clone = rownames(cloneStories)[i]
    xText = max(x) + ceiling(i/2)*1.5 - 0.9
    y0 = i/2 - floor(i/2)
    if ( labels ) {
      clone = names(cloneCols)[i]
      font = cloneLabels[[clone]]$font
      severity = cloneLabels[[clone]]$severity
      severity[font %in% c(2,4)] = -1   #put bold-fonted mutations at the top of the list
      ord = order(severity)
      muts = cloneLabels[[clone]]$label[ord]
      font = cloneLabels[[clone]]$font[ord]
      if ( length(muts) > 15 ) {
        muts = c(muts[1:15], 'and more...')
        font = c(font[1:15], 1)
      }
      col = cloneCols[i]
      text(rep(xText, length(muts)), y0 + 0.5 - (1:length(muts))/33, muts, col=col, adj=0, cex=0.9, font=font)
    }
  }

  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=5, col=rgb(0.7, 0.7, 0.7, 0.3))
  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=2, col=rgb(0.3, 0.3, 0.3, 0.3))
  text(1:ncol(stories), -0.02, colnames(stories), srt=20, cex=0.9)

  par(oma=rep(0, 4))
  par(mar=rep(4, 4))

  return(cloneCols)
}


#plots a set of parallel disjoint clones, and recurs to each clones subclones to be plotted on top.
addSubclone = function(cT, stories, ylims, dodgyness, colourPool=c(), margin=0.02, preNorm=1) {
  if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
  subClones = names(cT)
  subStories = stories[names(cT),,drop=F]
  maxSize = ylims[2,]-ylims[1,]
  cloneSum = margin*maxSize+colsums(t(margin*maxSize + t(subStories)))
  cloneSum[cloneSum == 0] = 1   #this happens if all subclones are 0 at a sample. this avoids NaNs.
  norm = pmin(1, maxSize/cloneSum)
  base = ylims[1,]
  usedCols = c()
  for ( subClone in subClones ) {
    subStory = subStories[subClone,]*norm
    range = rbind(base+margin*maxSize*norm, base + margin*maxSize*norm + subStory)
    range[2,] = pmax(base+margin*maxSize*norm, range[2,])
    subrange = rbind(range[1,], range[2]-margin*maxSize*norm)
    addStream(range, col=colourPool[1], dodgyness=dodgyness[subClone])
    usedCols = c(usedCols, colourPool[1])
    names(usedCols)[length(usedCols)] = subClone
    colourPool = colourPool[-1]
    if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
    if ( length(cT[[subClone]]) > 0 ) {
      out = addSubclone(cT[[subClone]], stories, range, dodgyness, colourPool, margin, preNorm=norm)
      usedCols = c(usedCols, out$usedCols)
      colourPool = out$colourPool
    }
    base = range[2,]
  }
  return(list('colourPool'=colourPool, usedCols=usedCols))
}

#plots the stream corresponding to a clone
addStream = function(ylims, col='grey', dodgyness=0) {
  for ( sample in 2:ncol(ylims) ) {
    x1 = sample-1
    x2 = sample
    y1h = ylims[2, sample-1]
    y1l = ylims[1, sample-1]
    y2h = ylims[2, sample]
    y2l = ylims[1, sample]
    range = c(0,1)
    if ( y1h == y1l & y2h > y2l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) < sample) )  {
      range[1] = 1-sqrt(y2h-y2l)
    }
    if ( y2h == y2l & y1h > y1l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) > sample)) {
      range[2] = sqrt(y1h-y1l)
    }
    addStreamSegment(x1, x2, y1l, y1h, y2l, y2h, range=range, col=col,dodgyness=dodgyness)
  }
}

#helper third degree polynomial
third = function(x, x0, a, b, c) a + b*(x-x0) + c*(x-x0)^3

#plots a smooth stream segment.
addStreamSegment = function(x1, x2, y1low, y1high, y2low, y2high, range=c(0,1), col, parts = 100, dodgyness=0) {
  sh = (y1high-y2high)/2
  sl = (y1low-y2low)/2
  z = (x1-x2)/2

  ah = (y1high+y2high)/2
  bh = 3*sh/(2*z)
  ch = -sh/(2*z*z*z)
  al = (y1low+y2low)/2
  bl = 3*sl/(2*z)
  cl = -sl/(2*z*z*z)
  x0 = (x1+x2)/2

  x = seq(from=x1, to=x2, length.out=parts)
  yhigh = third(x, x0, ah, bh, ch)
  ylow = third(x, x0, al, bl, cl)

  if ( range[1] == 0 & range[2] < 1 ) {
    xnorm = (x - x1)/(x2-x1)
    r = range[2]
    yShift = (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)^3*4
    yShift = ifelse(xnorm < r, yShift, 0)
    ylow = ylow + yShift
    yhigh = yhigh - yShift
    x = x[xnorm < r]
    ylow = ylow[xnorm < r]
    yhigh = yhigh[xnorm < r]
  }
  if ( range[1] > 0 & range[2] == 1 ) {
    xnorm = (x2 - x)/(x2-x1)
    r = (1-range[1])
    yShift = (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)^3*4
    yShift = ifelse(xnorm < r, yShift, 0)
    ylow = ylow + yShift
    yhigh = yhigh - yShift
    x = x[xnorm < r]
    ylow = ylow[xnorm < r]
    yhigh = yhigh[xnorm < r]
  }
  polygon(c(x, rev(x)), c(yhigh, rev(pmin(ylow, yhigh))), col=col, border=NA)
}

excludeClonesFromTree = function(tree, excludeClones) {
  if ( length(tree) == 0 ) return(tree)
  for ( clone in names(tree) ) {
    if ( clone %in% excludeClones )
      tree = c(tree[names(tree) != clone], excludeClonesFromTree(tree[[clone]], excludeClones))
    else
      tree[[clone]] = excludeClonesFromTree(tree[[clone]], excludeClones)
  }
  return(tree)
}

plotStories = function(stories, variants, col='default', lty='default', add=F, alpha=1, xlab='sample', ylab='clonality', lwd='default', errorBars=T, setPar=T, legend=T, labels=T, xlim='default', genome='hg19', xSpread=0.25, sampleOrder=NA,...) {
  if ( is.na(sampleOrder[1]) ) sampleOrder = colnames(stories$stories)
  clon = stories$stories[,sampleOrder,drop=F]
  ce = stories$errors[,sampleOrder,drop=F]

  Nsample = ncol(clon)
  Nmut = nrow(clon)

  if ( col[1] == 'default' | is.na(col[1]) | is.null(col[1]) | !(length(col) %in% c(1,Nmut)) ) {
    segcol = randomCols(row(clon)[,2:Nsample], a=alpha)
    errcol = randomCols(row(clon), a=alpha)
  }
  else if ( length(col) == 1 ) {
    col = rep(col, Nmut)
    segcol = rep(col, Nsample-1)
    errcol = rep(col, Nsample)
  }
  else if ( length(col) == Nmut ) {
    if ( all(rownames(clon) %in% names(col)))
      col = col[rownames(clon)]
    segcol = rep(col, Nsample-1)
    errcol = rep(col, Nsample)
  }

  if ( lty[1] == 'default' | !(length(lty) %in% c(1,Nmut)) ) {
    iMx = floor(1+(row(clon)-1)/8)
    seglty = randomLtys(iMx[,2:Nsample])
    errlty = randomLtys(iMx)
  }
  else if ( length(lty) == 1 ) {
    lty = rep(lty, Nmut)
    seglty = rep(lty, Nsample-1)
    errlty = rep(lty, Nsample)
  }
  else if ( length(lty) == Nmut ) {
    seglty = rep(lty, Nsample-1)
    errlty = rep(lty, Nsample)
  }

  
  if ( !add ) {
    if ( setPar ) {
      par(oma=rep(0, 4))
      par(mar=c(0, 4, 0, 0))
    }
    if ( xlim[1] == 'default' ) xlim = c(1, Nsample*1.3)
    plot(0,0, type='n', xlim=xlim, ylim=c(-.02,1), xlab=xlab, ylab=ylab, xaxt='n', frame.plot=F, ...)
  }
  if (lwd[1] == 'default' ) lwd = 1/(sqrt(0.1^2+ce[,1:(Nsample-1)]^2 + ce[,2:Nsample]^2)/0.2)^2
  segments(col(clon)[,1:(Nsample-1)]+(row(clon)[,1:(Nsample-1)] - 0.5 - Nmut/2)/Nmut*xSpread, clon[,1:(Nsample-1)],
           col(clon)[,2:Nsample    ]+(row(clon)[,2:Nsample    ] - 0.5 - Nmut/2)/Nmut*xSpread, clon[,2:Nsample    ],
           col=segcol, lwd=lwd, lty=seglty)
  if ( errorBars ) {
    if (lwd[1] == 'default' ) lwd = 1/(sqrt(0.1^2+ce^2)/0.2)^2
    segments(col(clon)+(row(clon) - 0.5 - Nmut/2)/Nmut*xSpread, noneg(clon - ce),
             col(clon)+(row(clon) - 0.5 - Nmut/2)/Nmut*xSpread, clon + ce,
             col=errcol, lwd=lwd, lty=errlty)
  }

  if ( !add ) {
    if ( legend ) {
      lbls = storyToLabel(stories, variants, genome)
      legCex = pmin(1, pmax(0.5, 45/nrow(lbls)))
      font = lbls$font
      severity = lbls$severity
      severity[font %in% c(2,4)] = -1   #put bold-fonted mutations at the top of the list
      ord = order(severity)
      legend('topright', lbls$label[ord], lwd=2, col=errcol[1:nrow(stories)][ord], lty=errlty[1:nrow(stories)][ord],
             cex= legCex, text.col=lbls$colour[ord], text.font=lbls$font[ord], seg.len=4)
    }
    if ( labels ) text(1:Nsample, -0.02, colnames(stories$stories), srt=20, cex=0.9)
    if ( setPar ) {
      par(oma=rep(0, 4))
      par(mar=rep(4, 4))
    }
  }
}

heatmapStories = function(stories, storyList, variants, col=NA, genome='hg19') {
  stories = stories[do.call(c, storyList),]
  if ( !is.na(col)[1] & all(names(storyList) %in% names(col)) ) {
    clone = do.call(c, lapply(1:length(storyList), function(i) rep(i, length(storyList[[i]]))))
    sideCol = do.call(c, lapply(1:length(storyList), function(i) rep(col[names(storyList)[i]], length(storyList[[i]]))))
  }
  else {
    clone = do.call(c, lapply(1:length(storyList), function(i) rep(i, length(storyList[[i]]))))
    sideCol = randomCols(clone)
  }
  
  clonalityMx = stories$stories

  labels = storyToLabel(stories, variants, genome)
  rownames(clonalityMx) = labels$label

  if ( nrow(clonalityMx) < 1000 ) {
    worked = try(makeHeatmap(clonalityMx, RowSideColors=sideCol, label='clonality'))
    if ( class(worked) == 'try-error' ) {
      catLog('Error in the heatmap. Trying without row/column clustering.\n')
      makeHeatmap(clonalityMx, RowSideColors=sideCol, Colv=NA, Rowv=NA)
    }
  }
  else {
    catLog('Too many stories for the built-in heatmap clustering, using default row ordering.\n')
    makeHeatmap(clonalityMx, RowSideColors=sideCol, Colv=NA, Rowv=NA)
  }
}

#' plots a heatmap
#'
#' @param mx numeric. The matrix to be plotted.
#' @param nCol numeric. The number of colours in the gradient. Default 200.
#' @param col character. The colour scheme to be used. The default 'default' gives a scale adapted for freqwuencies between 0 and 1, with extra attention to values close to 0. 'sunset' is a easier on the eye, but has less ability to resolve values close to the lowest value. 'DE' is meant for differential expression, and sets 0 to black, with red and blue gradient towards positive and negative values, saturating at the DEsaturation parameter.
#' @param maxVal numeric. Only for internal use. Default 'default'.
#' @param minVal numeric. Only for internal use. Default 'default'.
#' @param label character. The label for the colour scale on the side. Default ''.
#' @param DEsaturation numeric. The value where the colour saturates when col='DE'. Default log2(10).
#' @param ...  remaining arguments are passed to base heatmap(...)
#'
#' @details This function is a wrapped for the base heatmap() function. It got nicer default colours, doesn't normalise rows or columns by deafult, and has some support for differential expression data. Also prints a colour scale at the side.
#'
#'
#' @export
#' @examples
#' #random matrix to plot, centered around 0. Plot in 'DE' colours.
#' mx = matrix(rt(400, df=10), nrow=100)
#' makeHeatmap(mx, col='DE')
#' 
#' #random matrix to plot, between 0 and 1. Plot in default and sunset colours.
#' mx = matrix(runif(400), nrow=100)
#' makeHeatmap(mx)
#' makeHeatmap(mx, col='sunset')
#'
makeHeatmap = function(mx, nCol=200, col='default', maxVal='default', minVal='default', scale='none', label='', DEsaturation=log2(10), reverseGradient=F, ...) {
  if ( nrow(mx) < 2 | ncol(mx) < 2 ) {
    warning('Not two lines and two columns in the heatmap. skip.')
    return()
  }
  if ( maxVal == 'default' ) maxVal = max(mx, na.rm=T)
  if ( minVal == 'default' ) minVal = min(mx, na.rm=T)
  if ( col[1] == 'default' )
    col = colourGradient(cols=mcri(c('black', 'grey', 'cyan', 'blue', 'red')),
      anchors=c(0, 0.02, 0.2, 0.5, 1), steps=nCol, reverseGradient=reverseGradient)
  if ( col[1] == 'sunset' )
    col = colourGradient(cols=mcri(c('black', 'blue', 'cyan', 'orange')),
      anchors=c(0, 0.1, 0.5, 1), steps=nCol, reverseGradient=reverseGradient)
  DEmode = F
  if ( col[1] == 'DE' ) {
    DEmode = T
    zero = min(mx,na.rm=T)/(min(mx,na.rm=T)-max(mx,na.rm=T))
    upScale = min(1, DEsaturation/max(mx,na.rm=T))
    dnScale = min(1, -DEsaturation/min(mx,na.rm=T))
    if ( max(mx,na.rm=T) <= 0 )
      col = colourGradient(cols=mcri(c('cyan', 'blue', 'black')),
        anchors=c(max(0,min(1,zero - zero*dnScale)), max(0,min(1,zero - zero*dnScale/2)), max(0,min(1,zero))), steps=nCol,
        reverseGradient=reverseGradient)
    else if ( min(mx,na.rm=T) >= 0 )
      col = colourGradient(cols=mcri(c('black', 'red', 'orange')),
        anchors=c(max(0,min(1,zero)), max(0,min(1,zero + (1-zero)*upScale/2)), max(0,min(1,zero + (1-zero)*upScale))),
        steps=nCol, reverseGradient=reverseGradient)
    else
      col = colourGradient(cols=mcri(c('cyan', 'blue', 'black', 'red', 'orange')),
        anchors=c(max(0,min(1,zero - zero*dnScale)), max(0,min(1,zero - zero*dnScale/2)), max(0,min(1,zero)), max(0,min(1,zero + (1-zero)*upScale/2)), max(0,min(1,zero + (1-zero)*upScale))), steps=nCol, reverseGradient=reverseGradient)
    
  }

  ret = heatmap(mx, col=col, scale=scale, ...)
  
  barXmax = par('usr')[1]*0.88+par('usr')[2]*0.12
  barXmin = par('usr')[1]*0.92+par('usr')[2]*0.08
  barYmin = par('usr')[3]*0.8+par('usr')[4]*0.2
  barYmax = par('usr')[3]*0.2+par('usr')[4]*0.8
  lowY = barYmin + (0:(nCol-1) - 0.2)*(barYmax-barYmin)/nCol
  highY = barYmin + (1:nCol + 0.2)*(barYmax-barYmin)/nCol
  lowX = rep(barXmin, nCol)
  highX = rep(barXmax, nCol)
  barCols = col
  rect(lowX, lowY, highX, highY, col=barCols, border=NA)
  midTick = 0.5
  if ( DEmode ) midTick = -minVal/(maxVal-minVal)
  segments(rep(barXmin, 3), barYmin + c(0.002, midTick, 0.998)*(barYmax-barYmin),
           rep(barXmin-(barXmax-barXmin)*0.1, 3), barYmin + c(0.002, midTick, 0.998)*(barYmax-barYmin),
           lwd=2, col=barCols[c(1, round(nCol*midTick), nCol)])
  minDist = minVal
  maxDist = maxVal
  text(rep(barXmin-(barXmax-barXmin)*0.2, 3), barYmin + c(0.003, midTick, 0.997)*(barYmax-barYmin),
       c(round(c(minDist, minDist + (maxDist-minDist)*midTick, maxDist), 2)), adj=c(1, 0.5))
  text(barXmin-(barXmax-barXmin)*0.2, barYmin + 1.07*(barYmax-barYmin),
       label, adj=c(0.5, 0.5))

  invisible(ret)
}


#This function takes two colours (in a format that can be handled by col2rgb, such as "red", or rgb(0.1, 0.2, 0.3))
#and two weights and returns a weighted average between the two colours.
#Mainly a helper function for colourGradient, but can potentially find uses on itself.
combineColours = combineColors = function (col1, col2, w1=0.5, w2=0.5) {
  rgb1 = col2rgb(col1)
  rgb2 = col2rgb(col2)
  if ( w1 == Inf & w2 == Inf ) {
    w1 = 0.5
    w2 = 0.5
  }
  if ( w1 == Inf ) {
    w1 = 1
    w2 = 0
  }
  if ( w2 == Inf ) {
    w1 = 0
    w2 = 1
  }
  combinedRgb = (rgb1*w1 + rgb2*w2)/(w1+w2)
  combinedColour = do.call(rgb, as.list(pmin(1, pmax(0, combinedRgb/255))))
  return(combinedColour)
}
#Takes a vector of colours (in a format that can be handled by col2rgb, such as "red", or rgb(0.1, 0.2, 0.3))
#and an optional sorted vector of anchor points between 0 and 1. Returns a vector of colours of length @steps
#that gradually goes through the colours in the vector, hitting each colour at  fraction through the vector
#set by the anchor points. Defaults to a 100-step vector from blue through white to red.
colourGradient = colorGradient = function(cols=mcri(c('red', 'orange', 'white', 'cyan', 'blue')), steps=100, anchors=(1:length(cols)-1)/(length(cols)-1), reverseGradient=F ) {
  if ( length(anchors) != length(cols) ) {
    warning(paste0('colourGradient: Length of cols and anchors has to be the same. They are ', length(cols), ' and ', length(anchors), '. Returning default colour gradient.'))
  }
  N = length(cols)
  x = (0:(steps-1))/(steps-1)
  if ( any(anchors != sort(anchors)) ) warning('colourGradient: anchors not sorted. Sorting.')
  ord = order(anchors)
  anchors = anchors[ord]
  cols = cols[ord]
  anchors = c(-Inf, anchors, Inf)
  cols = c(cols[1], cols, cols[length(cols)])
  col2i = sapply(x, function(X) which(anchors > X)[1])
  col1i = col2i-1
  col1 = cols[col1i]
  col2 = cols[col2i]
  w1 = 1/abs(x-anchors[col1i])
  w2 = 1/abs(x-anchors[col2i])
  gradient = sapply(1:length(col1), function(i) combineColours(col1[i], col2[i], w1[i], w2[i]))
  if ( reverseGradient ) gradient = rev(gradient)
  names(gradient) = x
  return(gradient)
}


plotClones = function(cloneStories, allStories, storyList, variants, sampleOrder=NA, bg.al=0.1, maxError=1, ...) {
  if ( is.na(sampleOrder[1]) ) sampleOrder = colnames(cloneStories$stories)
  cloneStories$stories = cloneStories$stories[,sampleOrder]
  cloneStories$error = cloneStories$error[,sampleOrder]
  cloneNames = rownames(cloneStories$stories)

  cloneStories$stories = t(sapply(rownames(cloneStories), function(clone) {
    alls = allStories[storyList[[clone]],]
    return(colsums(alls$stories/alls$errors^2)/colsums(1/alls$errors^2))
    }))
  
  cloneCol = mcri(1:length(cloneNames)-1)
  fadedCloneCol = mcri(1:length(cloneNames)-1, al=bg.al)
  names(cloneCol) = names(fadedCloneCol) = cloneNames
  plotStories(cloneStories, variants, col=cloneCol, xSpread=0, errorBars=F, lwd=0, ...)  
  for ( clone in cloneNames ) {
    all = allStories[storyList[[clone]],]
    all = all[rowmeans(all$errors) < maxError,]
    cloneStories$stories = 
    plotStories(all, variants, add=T, col=fadedCloneCol[clone], lty=1, errorBars=F, xSpread=0)
  }
  plotStories(cloneStories, variants, col=cloneCol, xSpread=0, errorBars=F, lwd=5, add=T)  
}


reorderStories = function(stories, newOrder) {
  if ( !all(newOrder %in% colnames(stories$allStories$stories)) )
    stop("The new order ", do.call(paste, as.list(newOrder)), " does not use only sample names: ",
         do.call(paste, as.list(colnames(stories$allStories$stories))))
    
  stories$allConsistentStories$stories = stories$allConsistentStories$stories[,newOrder]
  stories$allConsistentStories$errors = stories$allConsistentStories$errors[,newOrder]
  stories$consistentClusteredStories$cloneStories$stories = stories$consistentClusteredStories$cloneStories$stories[,newOrder]
  stories$consistentClusteredStories$cloneStories$errors = stories$consistentClusteredStories$cloneStories$errors[,newOrder]
  stories$allStories$stories = stories$allStories$stories[,newOrder]
  stories$allStories$errors = stories$allStories$errors[,newOrder]
  stories$clusteredStories$stories = stories$clusteredStories$stories[,newOrder]
  stories$clusteredStories$errors = stories$clusteredStories$errors[,newOrder]
  return(stories)
}


cloneToCol = function(tree, colourPool=c()) {
  if ( length(colourPool) == 0 )
    colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange',
      'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
  
  subClones = names(tree)

  usedCols = c()
  for ( subClone in subClones ) {
    usedCols = c(usedCols, colourPool[1])
    names(usedCols)[length(usedCols)] = subClone
    colourPool = colourPool[-1]
    
    if ( length(colourPool) == 0 )
      colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange',
        'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
    if ( length(tree[[subClone]]) > 0 ) {
      out = cloneToCol(tree[[subClone]], colourPool)
      usedCols = c(usedCols, out$usedCols)
      colourPool = out$colourPool
    }
  }

  return(list('colourPool'=colourPool, usedCols=usedCols))
}
