
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
plotCNAheatmapWithGoI = function(GoI, Rdirectory, genome, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c()) {
  superFreq:::resetMargins()
  superFreq:::plotCNAbatchHeatmap(Rdirectory=Rdirectory, genome=genome,
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples)
  individuals = list.dirs(paste0(Rdirectory, "/"), full.names = F)
  individuals = individuals[individuals != ""]
  individuals = individuals[!(individuals %in% excludeIndividuals)]
  load(paste0(Rdirectory, '/', individuals[1], '/captureRegions.Rdata'))
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
plotCNAheatmapOverGoI = function(gene, Rdirectory, genome, padding=3e6, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c()) {
  individuals = list.dirs(paste0(Rdirectory, "/"), full.names = F)
  individuals = individuals[individuals != ""]
  individuals = individuals[!(individuals %in% excludeIndividuals)]
  load(paste0(Rdirectory, '/', individuals[1], '/captureRegions.Rdata'))

  crs = captureRegions[captureRegions$region == gene]
  maxx = chrToX(seqnames(crs[1]), max(end(crs)), genome=genome)
  minx = chrToX(seqnames(crs[1]), min(start(crs)), genome=genome)
  midx = (minx+maxx)/2

  superFreq:::resetMargins()
  superFreq:::plotCNAbatchHeatmap(Rdirectory=Rdirectory, genome=genome, xlim=c(minx-padding, maxx+padding),
                                  excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples)
  ymax = par("usr")[4]
  segments(minx, 0.5, minx, ymax-0.75, col=mcri('cyan'))
  segments(maxx, 0.5, maxx, ymax-0.75, col=mcri('cyan'))
  text(midx, 0.25, renameGoIs(gene, GoIakas), cex=0.8, col=mcri('cyan'))
}


#this is just a copy-paste from teh venetoclax AML batch. rewrite to be general.
plotZoomedCNAOnGoI = function(GoI, Rdirectory, individual, sample, genome, padding=3e6, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c()) {
  load(paste0(Rdirectory, '/', individual, '/clusters.Rdata'))
  
  
    cl = clusterList[[donor]]
    layout(matrix(1:2, ncol=1))
    cr = cl[[sample]]$CR
    use = xToChr(cr$x1, genome='hg38')=='20'
    cr = cr[use,]
    cluster = cl[[sample]]$clusters
    use = xToChr(cluster$x1, genome='hg38')=='20'
    cluster = cluster[use,]
    plot(xToPos((cr$x1+cr$x2)/2,genome='hg38'), cr$M, pch=16, cex=pmin(1,0.5/sqrt(cr$width)), ylim=c(-4,4), xlab='genomic position chr20', ylab='expression LFC compared to leucegene panel', main=sample)
    segments(-1e99, 0, 1e99, 0, col='grey')

    isCNN = cluster$call %in% c("AA", "AA?", "AA??")
    isAB = cluster$call == "AB"
    red = superFreq:::noneg(pmin(cluster$M, 1)) * (!isCNN) * 
      (!isAB)
    darkred = cluster$M > 1 | nchar(gsub("?", "", cluster$call)) > 
      4
    blue = superFreq:::noneg(pmin(-cluster$M, 1)) * (!isCNN) * 
      (!isAB)
    green = isCNN * (0.5 - cluster$f) * 2 * (!isAB)
    green = ifelse(is.na(green), 0, green)
    alpha = 0.1
    col = ifelse(darkred, mcri("red", 2*alpha), ifelse(red > 0, 
      mcri("red", alpha), ifelse(blue > 0, mcri("blue", alpha), 
                               mcri("green", alpha*(green>0)))))
    pos1 = xToPos(cluster$x1, genome='hg38')
    pos2 = xToPos(cluster$x2, genome='hg38')
    segments(pos1, cluster$M, pos2, cluster$M, lwd=2)
    posBCL = xToPos(midx, genome='hg38')
    segments(posBCL, -3, posBCL, 10, lwd=0.5)
    text(posBCL, -3.5, "BCL-XL", font=2, cex=0.8)
    rect(pos1, -10, pos2, 10, col=col, border=F)

    efreq = cl[[sample]]$eFreqs
    use = xToChr(efreq$x, genome='hg38')=='20'
    efreq = efreq[use,]
    plot(xToPos(efreq$x, genome='hg38'), efreq$var/efreq$cov, pch=16, cex=3/sqrt(efreq$cov), ylim=c(0,1), xlab='genomic position chr20', ylab='BAF')
    segments(-1e99, 0.5, 1e99, 0.5, col='grey')
    f = ifelse(cluster$postHet < 0.01, cluster$f, 0.5)
    segments(posBCL, 0.1, posBCL, 10, lwd=0.5)
    text(posBCL, 0.02, "BCL-XL", font=2, cex=0.8)
    segments(pos1, f, pos2, f, lwd=2)
    segments(pos1, 1-f, pos2, 1-f, lwd=2)
    rect(pos1, -10, pos2, 10, col=col, border=F)
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
plotCohortMutationHeatmap = function(GoI, Rdirectory, excludeIndividuals=c(), excludeSamples=c(), GoIakas=c()) {
  individuals = list.dirs(paste0(Rdirectory, "/"), full.names = F)
  individuals = individuals[individuals != ""]
  individuals = individuals[!(individuals %in% excludeIndividuals)]
  names(individuals) = individuals
  clusterList = lapply(individuals, function(ind) {
    load(paste0(Rdirectory, "/", ind, "/clusters.Rdata"))
    clusters = clusters[!(names(clusters) %in% excludeSamples)]
    return(clusters)
  })
  allVariantsList = lapply(individuals, function(ind) {
    load(paste0(Rdirectory, "/", ind, "/allVariants.Rdata"))
    allVariants$variants$variants = allVariants$variants$variants[!(names(allVariants$variants$variants) %in% excludeSamples)]
    return(allVariants)
  })
  sampleList = lapply(clusterList, names)
  
  #set up plot and sample/gene labels
  PID = setupCohortPlot(GoI, sampleList, GoIakas)

  #fill in CNAs as background colour for each box
  addCNAs(GoI, clusterList, sampleList)

  #add dots for point mutations
  addPointMutations(GoI, allVariantsList, sampleList)
}

#helper function
setupCohortPlot = function(GoI, sampleList, GoIakas=c()) {

  #set up scale and grid
  xScale = length(GoI)
  yScale = length(unlist(sampleList))
  xTicks = 1:xScale
  yTicks = 1:yScale

  #title
  main = 'placeholder title'
  
  #empty plot
  plot(0, type='n', xaxt='n', xlim=c(-0.15, 1)*(xScale+0.5), ylim=c(0, 1.3)*yScale, frame=F, yaxt='n', xlab='', ylab='', main=main)

  #grid
  segments(0:xScale + 0.5, 0.5, 0:xScale + 0.5, yScale*1.15, lwd=0.5, col='grey')
  segments(-0.03*xScale, 0:yScale + 0.5, xScale + 0.5, 0:yScale + 0.5, lwd=0.3, col=rgb(0.8, 0.8, 0.8))
  if ( length(sampleList) > 1 ) {
    indSeparators = cumsum(sapply(sampleList, length))
    indSeparators = indSeparators[-length(indSeparators)]
    segments(-0.15*xScale, indSeparators + 0.5, xScale + 0.5, indSeparators + 0.5, lwd=1, col='grey')
  }

  textCol = rgb(0.5,0.5,0.5)
  for ( ind in names(sampleList) ) {
    sys = cumsum(sapply(sampleList, length))[ind] - length(sampleList[[ind]]) + 1:length(sampleList[[ind]])
    text(-0.03*xScale, sys, sampleList[[ind]], col=textCol, adj=c(1, 0.5), cex=0.8, font=2)
  }

  #topbar
  text(xTicks, yScale+1, renameGoIs(GoI, GoIakas), adj=c(0, 0.5), srt=90, cex=0.8, font=2, col=textCol)

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
        col = callToColHeatmap(call, sqrt(clonality))
        x = xs[gene]
        rect(x-0.45, y-0.5, x+0.45, y+0.5, border=F, col=col)
      }
    }
  }
  
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
addPointMutations = function(GoI, allVariantsList, sampleList) {
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
      q = allVariantsList[[ind]]$variants$variants[[sample]]
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
        bordercol = 'white'
        if ( any(qg$severity < 10) ) bordercol = mcri('orange')
        x = xs[gene]
        cex=sqrt(cex*1.2)
        points(x, y, col='white', cex=cex, pch=ifelse(qg$severity < 10, 24, 21), lwd=1.5, bg='black')
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


findGoI = function(Rdirectory='R', cosmicCensus=T, clinvarPathogenic=T, clinvarAny=F, excludeIndividuals=c(), excludeSamples=c(), maxNumberOfGenes=50) {
  qsList = superFreq:::loadQsList(Rdirectory=Rdirectory, excludeIndividuals=excludeIndividuals, excludeSamples=excludeSamples)
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
      (isCosmicCensus*cosmicCensus | isClinvarPathogenic*clinvarPathogenic | isClinvarAny*clinvarAny)
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
