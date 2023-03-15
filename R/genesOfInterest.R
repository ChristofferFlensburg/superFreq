
#' plots a heatmap of copy numbers, highlighting genes of interest
#'
#' @param GoI character vector. Genes of interest.
#' @param Rdirectory character The Rdirectory where the superFreq run was run.
#' @param genome character. The genome the sample is aligned to. 'hg19', 'hg38' or 'mm10'
#'
#' @details This function plots a heatmap of copy numbers across the samples and genome. Genes of Interest are highlighted. Can be run after the superFreq run is completed.
#'
#'
#' @export
#'
plotCNAheatmapWithGoI = function(GoI, Rdirectory, genome, metaDataFile, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c(), cpus=1) {
  superFreq:::resetMargins()
  superFreq:::plotCNAbatchHeatmap(Rdirectory=Rdirectory, metaDataFile=metaDataFile, genome=genome,
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  captureRegions = loadAnyCaptureRegions(Rdirectory, metaDataFile, excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples)
  ymax = par("usr")[4]
  for ( gene in GoI ) {
    crs = captureRegions[captureRegions$region == gene]
    maxx = chrToX(seqnames(crs[1]), max(end(crs)), genome=genome)
    minx = chrToX(seqnames(crs[1]), min(start(crs)), genome=genome)
    midx = (minx+maxx)/2
    segments(midx, 0.5, minx, ymax-0.75, col=mcri('cyan', 0.5), lwd=0.5)
    text(midx, 0.25, renameGoIs(gene, GoIakas), cex=0.8, col=mcri('cyan'))
  }
}

#' plots a heatmap of copy numbers, focusing on a gene of interest
#'
#' @param GoI character vector. Genes of interest.
#' @param Rdirectory character The Rdirectory where the superFreq run was run.
#' @param genome character. The genome the sample is aligned to. 'hg19', 'hg38' or 'mm10'
#' @param padding numeric. How many basepairs next to the gene to show. default 3e6.
#'
#' @details This function plots a heatmap of copy numbers across the samples, focused at a gene of interest. Can be run after the superFreq run is completed.
#'
#'
#' @export
#'
plotCNAheatmapOverGoI = function(gene, Rdirectory, genome, metaDataFile, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c(), cpus=1, padding=3e6) {
  captureRegions = loadAnyCaptureRegions(Rdirectory, metaDataFile, excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples)
  crs = captureRegions[captureRegions$region == gene]
  maxx = chrToX(seqnames(crs[1]), max(end(crs)), genome=genome)
  minx = chrToX(seqnames(crs[1]), min(start(crs)), genome=genome)
  midx = (minx+maxx)/2
  
  superFreq:::resetMargins()
  superFreq:::plotCNAbatchHeatmap(Rdirectory=Rdirectory, metaDataFile=metaDataFile, genome=genome, xlim=c(minx-padding, maxx+padding),
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  ymax = par("usr")[4]
  segments(minx, 0.5, minx, ymax-0.75, col=mcri('cyan'))
  segments(maxx, 0.5, maxx, ymax-0.75, col=mcri('cyan'))
  text(midx, 0.25, renameGoIs(gene, GoIakas), cex=0.8, col=mcri('cyan'))
}



#' plots mutations in samples in Genes of Interest
#'
#' @param GoI character vector. Genes of interest for the analysis
#' @param Rdirectory character The Rdirectory where the superFreq run was run.
#'
#' @details This function plots a mutation matrix showing CNAs and point mutations in the Genes of Interest across the analysed samples in the Rdirectory. Should be run after the superFreq() run has completed.
#'
#' @export
#'
plotCohortMutationHeatmap = function(GoI, Rdirectory, metaDataFile, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c(), cpus=1) {
  clusterList = superFreq:::loadClusterList(Rdirectory=Rdirectory, metaDataFile=metaDataFile, excludeIndividuals=excludeIndividuals,
                                excludeSamples=excludeSamples, cpus=cpus)
  sampleList = lapply(clusterList, names)
  qsList = superFreq:::loadQsList(Rdirectory=Rdirectory, metaDataFile=metaDataFile,
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  
  #set up plot and sample/gene labels
  superFreq:::setupCohortPlot(GoI, sampleList, GoIakas)

  #fill in CNAs as background colour for each box
  colMx = superFreq:::addCNAs(GoI, clusterList, sampleList)

  #add dots for point mutations
  superFreq:::addPointMutations(GoI, qsList, sampleList, colMx=colMx)
}

#helper function
setupCohortPlot = function(GoI, sampleList, GoIakas=c()) {

  #set up scale and grid
  xScale = length(GoI)
  yScale = length(unlist(sampleList))
  xTicks = 1:xScale
  yTicks = 1:yScale

  #title
  main = 'mutation matrix'
  
  #empty plot
  plot(0, type='n', xaxt='n', xlim=c(-0.15, 1)*(xScale+0.5), ylim=c(0, 1.3)*yScale, frame=F, yaxt='n', xlab='', ylab='', main=main)

  #grid
  segments(0:xScale + 0.5, 0.5, 0:xScale + 0.5, yScale*1.15, lwd=0.5, col='grey')
  if ( yScale < 100 )
    segments(-0.03*xScale, 0:yScale + 0.5, xScale + 0.5, 0:yScale + 0.5, lwd=0.3*min(1, 30/yScale), col=rgb(0.8, 0.8, 0.8))
  if ( length(sampleList) > 1 ) {
    indSeparators = cumsum(sapply(sampleList, length))
    indSeparators = indSeparators[-length(indSeparators)]
        #if many samples, thinner separating lines or no lines at all.
    if ( length(sampleList) < 100 )
      segments(-0.15*xScale, indSeparators + 0.5, xScale + 0.5, indSeparators + 0.5, lwd=min(1, 30/length(sampleList)), col='grey')
  }

  textCol = rgb(0.5,0.5,0.5)
  for ( ind in names(sampleList) ) {
    sys = cumsum(sapply(sampleList, length))[ind] - length(sampleList[[ind]]) + 1:length(sampleList[[ind]])
    text(-0.03*xScale, sys, sampleList[[ind]], col=textCol, adj=c(1, 0.5), cex=0.8, font=2)
  }

  #topbar
  text(xTicks, yScale+1, superFreq:::renameGoIs(GoI, GoIakas), adj=c(0, 0.5), srt=90, cex=0.8, font=2, col=textCol)

  #legend
  legend('topleft', c('ampli', 'gain', 'LOH', 'loss', 'SNV', 'biall', 'trunc'), bg='white', pch=c(15, 15, 15, 15, 16, 16, 17), col=mcri(c('darkred', 'red', 'green', 'blue', 'black', 'black', 'black')), pt.cex=c(2,2,2,2,1,1.4,1.2), pt.lwd=c(1,1,1,1,1,1,2))
}

#helper function
addCNAs = function(GoI, clusterList, sampleList) {
  #set up scale and grid
  xScale = length(GoI)
  yScale = length(unlist(sampleList))
  xTicks = 1:xScale
  yTicks = 1:yScale

  cr = clusterList[[1]][[1]]$CR
  XoI = (cr[GoI,]$x1+cr[GoI,]$x2)/2
  names(XoI) = GoI
  xs = 1:length(GoI)
  names(xs) = GoI

  colMx = matrix('black', ncol=length(GoI), nrow=length(unlist(sampleList)))
  rownames(colMx) = unlist(sampleList)
  colnames(colMx) = GoI
  for ( ind in names(sampleList) ) {
    sys = cumsum(sapply(sampleList, length))[ind] - length(sampleList[[ind]]) + 1:length(sampleList[[ind]])
    samples = sampleList[[ind]]
    names(sys) = samples
    for ( sample in samples ) {
      clusters = clusterList[[ind]][[sample]]$clusters
      y = sys[sample]
      for ( gene in GoI ) {
        call = clusters$call[clusters$x2 > XoI[gene] & clusters$x1 < XoI[gene]]
        if  ( length(call) == 0 | call %in% c('AB', 'AB?', 'AB??') ) next
        clonality = clusters$clonality[clusters$x2 > XoI[gene] & clusters$x1 < XoI[gene]]
        col = superFreq:::callToColHeatmap(call, sqrt(clonality))
        colMx[sample, gene] = col
        x = xs[gene]
        rect(x-0.45, y-0.5, x+0.45, y+0.5, border=F, col=col)
      }
    }
  }

  return(invisible(colMx))
}

#helper function
callToColHeatmap = function(call, clonality) {
  if ( call %in% c('A', 'A?', 'A??', 'B', 'B?', 'B??') ) return(mcri('blue', clonality))
  if ( call %in% c('AAB', 'AAB?', 'AAB??', 'BBA', 'BBA?', 'BBA??') ) return(mcri('red', clonality))
  if ( call %in% c('AA', 'AA?', 'AA??', 'BB', 'BB?', 'BB??', 'AAA', 'AAA?', 'AAA??') ) return(mcri('green', clonality))
  if ( call %in% c('CL', 'CL?', 'CL??') ) return(mcri('darkblue', clonality))
  if ( nchar(gsub('?', '', call)) > 3 ) {
    effectiveN = nchar(gsub('?', '', call))*clonality
    return(mcri('darkred', pmin(1, effectiveN/5)))
  }
  stop('Couldnt find colour.')
}


#helper function
addPointMutations = function(GoI, qsList, sampleList, colMx) {
  #set up scale and grid
  xScale = length(GoI)
  yScale = length(unlist(sampleList))
  xTicks = 1:xScale
  yTicks = 1:yScale

  xs = 1:length(GoI)
  names(xs) = GoI

  for ( ind in names(sampleList) ) {
    sys = cumsum(sapply(sampleList, length))[ind] - length(sampleList[[ind]]) + 1:length(sampleList[[ind]])
    samples = sampleList[[ind]]
    names(sys) = samples
    for ( sample in samples ) {
      q = qsList[[ind]][[sample]]
      q = q[q$somaticP > 0.5 & q$severity < 11 & q$inGene %in% GoI & q$var > 0.15*q$cov,]
      y = sys[sample]
      for ( gene in GoI ) {
        qg = q[q$inGene == gene,]
        if  ( nrow(qg) == 0 ) next
        cex=1
        if ( nrow(qg) > 1 ) cex=2
        else if ( qg$var > qg$cov*0.6 & pbinom(qg$var, qg$cov, 0.5, lower.tail=F) < 0.05 ) cex=2
        else if ( qg$var < qg$cov/4 ) cex = 4*qg$var/qg$cov
        pointcol = 'black'
        if ( yScale >= 100 ) pointcol = colMx[sample, gene]
        bordercol = 'white'
        x = xs[gene]
        cex=sqrt(cex*1.2)*min(1, max(0.6, sqrt(50/yScale)))
        lwd = 1.5*min(1, max(0.6, sqrt(50/yScale)))
        points(x, y, col=bordercol, cex=cex, pch=ifelse(qg$severity < 10, 24, 21), lwd=lwd, bg=pointcol)
      }
    }
  }
}


renameGoIs = function(genes, GoIakas) {
  for ( SYMBOLname in names(GoIakas) ) {
    genes = gsub(SYMBOLname, GoIakas[SYMBOLname], genes)
  }
  return(genes)
}


findGoI = function(Rdirectory, metaDataFile, cosmicCensus=T, clinvarPathogenic=T, clinvarAny=F, excludeIndividuals=c(), excludeSamples=c(), maxNumberOfGenes=50, cpus=1, genome='hg19') {
  qsList = superFreq:::loadQsList(Rdirectory=Rdirectory, metaDataFile=metaDataFile, excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  qs = do.call(c, qsList)
  GoIlist = lapply(qs, function(q) {
    isSomatic = q$somaticP > 0.5
    isCoding = q$severity < 11
    isntSubclonal = q$var > 0.15*q$cov
    isCosmicCensus = q$isCosmicCensus
    isClinvarAny = q$ClinVar_ClinicalSignificance != ''
    isClinvarPathogenic =
      grepl('[pP]athogenic,', paste0(q$ClinVar_ClinicalSignificance)) |
      grepl('[pP]athogenic$', paste0(q$ClinVar_ClinicalSignificance))
    use =
      isSomatic & isCoding & isntSubclonal &
      ((isCosmicCensus&cosmicCensus) | (isClinvarPathogenic&clinvarPathogenic) | (isClinvarAny&clinvarAny))
    if ( genome == 'mm10' )
    	use = isSomatic & isCoding & isntSubclonal & (isCosmicCensus | !cosmicCensus)

    return(unique(q$inGene[use]))
    })
  GoI = unique(unlist(GoIlist))
  hits = sapply(GoI, function(gene) {
    sum(sapply(names(qsList), function(ind) {
      qs = qsList[[ind]]
      GoIsublist = GoIlist[paste0(ind, '.', names(qs))]
      genesInInd = unique(unlist(GoIsublist))
      hit = gene %in% genesInInd
      return(hit)
    }))
  })
  GoI = GoI[order(-hits)]
  if ( length(GoI) > maxNumberOfGenes ) GoI = GoI[1:maxNumberOfGenes]
}




#' returns relevant mutations in samples in Genes of Interest
#'
#' @param GoI character vector. Genes of interest for the analysis
#' @param Rdirectory character The Rdirectory where the superFreq run was run.
#'
#' @details This function plots a mutation matrix showing CNAs and point mutations in the Genes of Interest across the analysed samples in the Rdirectory. Should be run after the superFreq() run has completed.
#'
#' @export
#'
getCohortMutationMatrix = function(GoI, Rdirectory, metaDataFile, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c(), cpus=1, includeGermlineLike=T) {
  clusterList = superFreq:::loadClusterList(Rdirectory=Rdirectory, metaDataFile=metaDataFile, excludeIndividuals=excludeIndividuals,
                                excludeSamples=excludeSamples, cpus=cpus)
  sampleList = lapply(clusterList, names)
  qsList = superFreq:::loadQsList(Rdirectory=Rdirectory, metaDataFile=metaDataFile,
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples, cpus=cpus)
  

  #get CNA calls and clonality
  cnaMx = superFreq:::getCohortCNAmatrix(GoI, clusterList, sampleList)

  #add dots for point mutations
  snvMx = getCohortSNVmatrix(GoI, qsList, sampleList, colMx=colMx, includeGermlineLike=includeGermlineLike)

  return(list(cnaMx=cnaMx, snvMx=snvMx))
}

#helper function
getCohortSNVmatrix = function(GoI, qsList, sampleList, colMx, includeGermlineLike=includeGermlineLike) {

  biallelicMx = matrix(FALSE, ncol=length(unlist(sampleList)), nrow=length(GoI))
  mostSevereMx = matrix('', ncol=length(unlist(sampleList)), nrow=length(GoI))
  rownames(biallelicMx) = rownames(mostSevereMx) = GoI
  colnames(biallelicMx) = colnames(mostSevereMx) = unlist(sampleList)
  
  for ( ind in names(sampleList) ) {
    samples = sampleList[[ind]]
    for ( sample in samples ) {
      q = qsList[[ind]][[sample]]
      q = q[q$somaticP > 0.5 & q$severity < 11 & q$inGene %in% GoI & q$var > 0.15*q$cov,]
	  if ( !includeGermlineLike ) q = q[!q$germline | is.na(q$germline),]
      for ( gene in GoI ) {
        qg = q[q$inGene == gene,]
        if  ( nrow(qg) == 0 ) next
        if ( nrow(qg) > 1 ) biallelicMx[gene, sample] = TRUE
        else biallelicMx[gene, sample] = qg$var > qg$cov*0.6 & pbinom(qg$var, qg$cov, 0.5, lower.tail=F) < 0.05
        mostSevereMx[gene, sample] = qg$type[which.min(qg$severity)]
      }
    }
  }

  return(list(biallelicMx=biallelicMx, mostSevereMx=mostSevereMx))  
}

#helper function
getCohortCNAmatrix = function(GoI, clusterList, sampleList) {
  #set up scale and grid

  cr = clusterList[[1]][[1]]$CR
  XoI = (cr[GoI,]$x1+cr[GoI,]$x2)/2
  names(XoI) = GoI

  callMx = matrix('', ncol=length(unlist(sampleList)), nrow=length(GoI))
  clonalityMx = matrix(0, ncol=length(unlist(sampleList)), nrow=length(GoI))
  rownames(callMx) = rownames(clonalityMx) = GoI
  colnames(callMx) = colnames(clonalityMx) = unlist(sampleList)
  for ( ind in names(sampleList) ) {
    samples = sampleList[[ind]]
    for ( sample in samples ) {
      clusters = clusterList[[ind]][[sample]]$clusters
      for ( gene in GoI ) {
        call = clusters$call[clusters$x2 > XoI[gene] & clusters$x1 < XoI[gene]]
        if  ( length(call) == 0 ) next
        clonalityMx[gene, sample] = clusters$clonality[clusters$x2 > XoI[gene] & clusters$x1 < XoI[gene]]
        callMx[gene, sample] = call
      }
    }
  }

  return(list(callMx=callMx, clonalityMx=clonalityMx))
}
