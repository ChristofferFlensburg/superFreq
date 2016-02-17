

#This function outputs jpg plots to the target directory summarising somatic mutations
#and CNVs over the genome for all samples
makeSummaryPlot = function(variants, cnvs, normals, individuals, timePoints, plotDirectory, genome='hg19', cpus=1, forceRedo=F, plotName='summary') {
  summaryDirectory = paste0(plotDirectory, '/summary/')
  if ( !file.exists(summaryDirectory) ) dir.create(summaryDirectory)
  plotFile = paste0(summaryDirectory, '/', plotName, '.png')
  if ( !file.exists(plotFile) | forceRedo ) {
    catLog('Making summary plot..')
    png(plotFile, width=20, height=10, res=300, units='in')
    plotSummary(variants, cnvs, normals, individuals, timePoints, genome, cpus)
    dev.off()
    for ( chr in names(chrLengths(genome)) ) {
      catLog(chr, '..', sep='')
      xmin = cumsum(chrLengths(genome))[chr]-chrLengths(genome)[chr]
      xmax = cumsum(chrLengths(genome))[chr]
      plotFile = paste0(summaryDirectory, '/', plotName, chr, '.jpg')
      jpeg(plotFile, width=20, height=10, res=300, units='in')
      plotSummary(variants, cnvs, normals, individuals, timePoints, genome, cpus, xmin=xmin, xmax=xmax)
      dev.off()
    }
    catLog('done.\n')
  }
  else catLog('Summary plot already in place, not redoing.\n')
}


#helper function that does all the work.
plotSummary = function(variants, cnvs, normals, individuals, timePoints, genome, cpus, xmin=0, xmax=sum(chrLengths(genome))) {
  N = length(cnvs)
  variants$variants = lapply(variants$variants, function(q) q[q$x > xmin & q$x < xmax,])
  cnvs = lapply(cnvs, function(cnv) {
    list(
         'clusters'=cnv$clusters[cnv$clusters$x2 > xmin & cnv$clusters$x1 < xmax,],
         'CR'=cnv$CR[cnv$CR$x2 > xmin & cnv$CR$x1 < xmax,],
         'freqs'=cnv$CR[cnv$freqs$x > xmin & cnv$freqs$x < xmax,])
  })
  if ( all(lapply(cnvs, function(cnv) nrow(cnv$CR)) == 0) ) return()
  
  #set up plot
  iToY = rank(individuals, ties.method='first')
  yToI = sapply(1:N, function(y) which(iToY == y))
  par(oma=rep(0, 4))
  par(mar=rep(0, 4))
  plot(0, type='n', xlim=c(xmin-(xmax-xmin)*0.02, xmax + (xmax-xmin)*0.2), ylim = c(-1.5,N+0.5), xaxt='n', yaxt='n', frame.plot=F,
       xlab='', ylab='')
  segments(c(rep(xmin, N+1), xmin, xmax), c(0:N, 0, 0), c(rep(xmax, N+1), xmin, xmax), c(0:N, N, N),
           lwd=c(5, ifelse(individuals[yToI]==c(individuals[yToI[2:N]], ''), 0.5, 4), 5, 5), col='grey')
  addChromosomeLines(ylim=c(0, N*1.02), col='grey', lwd=1.4, genome=genome)
  text(xmin-(xmax-xmin)*0.03, 1:N-0.35, names(cnvs)[yToI], cex=0.7)

  #plot the CNVs
  for (n in 1:N ) {
    cnv = cnvs[[n]]$clusters[!(cnvs[[n]]$clusters$call %in% c('AB', 'AB?', 'AB??')),]
    if ( nrow(cnv) == 0 ) next
    col = rgb(pmin(1, pmax(0, 1-cnv$M)), pmin(1, pmax(0, 1-abs(cnv$M))), pmin(1, pmax(0, 1+cnv$M)))
    rect(cnv$x1, iToY[n]-0.3 - pmin(1, abs(cnv$M))*0.2,
         cnv$x2, iToY[n]-0.3 + pmin(1, abs(cnv$M))*0.2, col=col, border=NA)
    cnv = cnvs[[n]]$clusters[cnvs[[n]]$clusters$call == 'AA',]
    if ( nrow(cnv) == 0 ) next
    col = rgb(1-cnv$clonality, 1-0.3*(cnv$clonality), 1-cnv$clonality)
    rect(cnv$x1, iToY[n]-0.3 - cnv$clonality*0.2,
         cnv$x2, iToY[n]-0.3 + cnv$clonality*0.2, col=col, border=NA)
  }

  #plot point mutations
  genes = list()
  #Loop over individuals
  for (ind in unique(individuals) ) {
    if ( sum(individuals == ind) == 1 ) {
      n = which(individuals == ind)[1]
      q = variants$variants[[n]]
      use = q$somaticP > 0.5
      if ( 'severity' %in% names(q) )
        use = q$somaticP > 0.5 &
             (is.na(q$germline) | !q$germline) &
             (is.na(q$severity) | q$severity < 11)
      q = q[use,]
      points(q$x, rep(iToY[n]-0.75, sum(use)), pch=19, cex=1.2*sqrt(q$var/q$cov), col=rgb(0,0,0,pmin(1,sqrt(q$cov/100))))
      genes[[names(cnvs)[n]]] = unique(variants$SNPs[as.character(q$x),]$inGene)
      next
    }
    #set up information common for all indiviuals
    qs = variants$variants[individuals == ind]
    var = do.call(cbind, lapply(qs, function(q) q$var))
    keep = rowSums(var) > 0
    qs = lapply(qs, function(q) q[keep,])
    fs = do.call(cbind, lapply(qs, function(q) q$var/q$cov))
    cov = do.call(cbind, lapply(qs, function(q) q$cov))
    var = do.call(cbind, lapply(qs, function(q) q$var))
    flag = do.call(cbind, lapply(qs, function(q) q$flag))
    f = rowSums(var)/rowSums(cov)
    fs[is.na(fs)] = -0.02

    #loop over samples for each individual
    for ( col in 1:ncol(var) ) {
      use = qs[[col]]$somaticP > 0.5 & (!qs[[col]]$germline | is.na(qs[[col]]$germline))
      if ( sum(use) == 0 ) next
      if ( 'severity' %in% names(qs[[col]]) )
        use = qs[[col]]$somaticP > 0.5 &
             (is.na(qs[[col]]$germline) | !qs[[col]]$germline) &
             (is.na(qs[[col]]$severity) | qs[[col]]$severity < 11)
      
      p = unlist(mclapply(which(use), function(row)
        pBinom(cov[row,col], var[row,col], sum(var[row,-col])/sum(cov[row,-col])),
        mc.cores=cpus))
      significant = p < 1/pmax(20, sum(use)^0.75)
      q = qs[[col]][use,]
      uniqueVar = rowsums(var[use,-col]) == 0
      smallerVar = (rowsums(var[use,-col])/rowsums(cov[use,-col])) > fs[use,col]
      cols = ifelse(significant, ifelse(uniqueVar, 'red', ifelse(smallerVar, 'blue', 'orange')) ,rgb(0,0,0,pmin(1,sqrt(q$cov/100))))
    
      #colour SNPs that are different between individuals
      #red for unique, orange for increasing freq, blue for decreasing freq.
      ord = order(significant + uniqueVar + smallerVar/2)
      n = which(individuals == ind)[col]
      points(q$x[ord], rep(iToY[n]-0.75, sum(use))[ord], pch=19, cex=1.2*sqrt(q$var/q$cov)[ord], col=cols[ord])

      genes[[names(cnvs)[n]]] = unique(variants$SNPs[as.character(q$x),]$inGene)
    }
  }

  #Highlight reccuring genes
  allGenes = Reduce(union, genes)
  geneCounts = sapply(allGenes, function(gene) length(unique(individuals[sapply(genes, function(gs) gene %in% gs)])))
  names(geneCounts) = allGenes
  i = 0
  print = allGenes[geneCounts > i]
  while( length(print) > 50 & sum(geneCounts > i+1) > 5 ) {
    i = i+1
    print = allGenes[geneCounts > i]
  }
  if ( length(print) > 0 ) {
    for ( name in names(cnvs) ) {
      q = variants$variants[[name]]
      
      use = q$somaticP > 0.5 & (!q$germline | is.na(q$germline))
      if ( 'severity' %in% names(q) )
        use = q$somaticP > 0.5 &
             (is.na(q$germline) | !q$germline) &
             (is.na(q$severity) | q$severity < 11)

      is = which(variants$SNPs$x %in% q$x) 
      use[use] = variants$SNPs[is,][as.character(q$x[use]),]$inGene %in% print
      if ( sum(use) == 0 ) next
      q = q[use,]
      n = which(names(cnvs) == name)
      segments(q$x, rep(iToY[n]-0.6, sum(use)), q$x, rep(iToY[n]-0.9, sum(use)), lwd=0.5)
    }
    
    geneX = sapply(print, function(gene) mean(variants$SNPs$x[variants$SNPs$inGene == gene]))
    x = sort(geneX)
    #textX = seq(from=xmin, to=xmax, along.with=print)
    textX = -spreadPositions(-x, (xmax-xmin)/50)
    textY = rep(c(-0.55, -0.85, -1.15, -1.45), length(print))[1:length(print)] - 0.2
    segments(textX, textY+0.2, x, -0.1, col=rgb(0.8, 0.8, 0.8))
    label = paste0(print, ' (', geneCounts[print], ')')
    text(textX, textY, label[order(geneX)], cex=0.7)
    text(xmax + (xmax-xmin)*0.12, -1, paste0('<-- somatic SNV in more than ',i,' patients.'))
  }

  legend('right', c('somatic point mutation', '50% frequency', '100% frequency', 'new mutation', 'increasing mutation', 'decreasing mutation',
                    'copy number gain', 'small gain', 'copy number loss', 'small loss', 'CNN LoH', 'small CNN LoH', 'recurring mutated gene'),
         pch=c(19, 19, 19, 19, 19, 19, 15, 15, 15, 15, 15, 15, 108),
         col=c(rgb(0.8, 0.8, 0.8), rgb(0.8, 0.8, 0.8), rgb(0.8, 0.8, 0.8), 'red', 'orange', 'blue', 'blue', rgb(0.7, 0.7, 1), 'red', rgb(1, 0.7, 0.7), rgb(0, 0.7, 0), rgb(0.7, 0.91, 0.7), 'black'),
         pt.cex=c(1,1.2*sqrt(0.5),1.2,1,1,1, 2, 1, 2, 1, 2, 1, 1), bg='white')
  
  if ( length(print) > 0 ) invisible(print[order(geneX)])
  invisible(c())
}

#spreads out positions. useful for plotting labels
spreadPositions = function(x, d) {
  ord = order(x)
  prevX = -Inf
  for ( i in ord ) {
    if ( x[i] < prevX+d ) x[i] = prevX+d
    prevX = x[i]
  }
  return(x)
}

#scalar product of the two normalsied vectors x and y.
scalarNorm = function(x, y) {
  return(sum(x*y)/norm(x)/norm(y))
}
