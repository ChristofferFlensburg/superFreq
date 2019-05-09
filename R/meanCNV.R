projectMeanCNV = function(metaData, project, cpus=1, cosmicDirectory='', onlyDNA=T, clonalityCut=0.4, includeNormal=F, cnvWeight=1, forceRedoMean=F, forceRedoVariants=F, forceRedoMeanPlot=F, forceRedoMatrixPlot=F, genome='hg19', ignoreCNAonly=F) {
  samples = inProject(metaData, project, includeNormal=includeNormal, onlyDNA=onlyDNA)
  if  ( length(samples) == 0 ) {
    warning('No cancer to analyse in project', project)
    return()
  }
  individuals = metaData$samples[samples,]$INDIVIDUAL

  variants = getProjectVariants(metaData, project, cpus=cpus, onlyDNA=onlyDNA,
    includeNormal=includeNormal, forceRedo=forceRedoVariants)
  
  saveFile = paste0(metaData$project[project,]$Rdirectory, '/meanCNV.Rdata')
  if ( file.exists(saveFile) & !forceRedoMean ) {
    catLog('loading cohort data...')
    load(saveFile)
    catLog('done.\n')
  }
  else {
    meanCNV = getMeanCNV(metaData, samples, variants, genome=genome, clonalityCut=clonalityCut, cpus=cpus, ignoreCNAonly=ignoreCNAonly)
    
    catLog('Saving meanCNV...')
    save(meanCNV, file=saveFile)
    catLog('done.\n')
  }

  catLog('Plotting mean CNV...')
  plotMeanCNVtoFile(metaData, project, meanCNV, cosmicDirectory=cosmicDirectory, forceRedo=forceRedoMeanPlot, genome=genome)
  catLog('mutation matrix...')
  plotMultipageMutationMatrix(metaData, meanCNV, project, cosmicDirectory=cosmicDirectory, nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot, genome=genome)
  catLog('done.\n')

  variants$variants = lapply(variants$variants, function(q) q[q$severity < 11,])
  outputSomaticVariants(variants, genome=genome, plotDirectory=metaData$project[project,]$plotDirectory, cpus=cpus)
  
  subgroups = getSubgroups(metaData, project, includeNormal=includeNormal)
  catLog('Subgroups:\n')
  for ( sg in subgroups ) cat(sg, '\n')
  if ( length(subgroups) > 0 ) { 
    for ( subgroup in subgroups ) {
      catLog('Plotting mean CNV for subgroup', subgroup, '...')
      plotDirectory = paste0(metaData$project[project,]$plotDirectory, '/', subgroup)
      meanCNVs = plotSubgroupCNVtoFile(metaData, meanCNV, project, subgroup, includeNormal=includeNormal, cosmicDirectory=cosmicDirectory,
        forceRedo=forceRedoMeanPlot, genome=genome)
      plotMeanCNVtoFile(metaData, project, meanCNVs$inGroup, cosmicDirectory=cosmicDirectory, forceRedo=forceRedoMeanPlot,
                        plotFile=paste0(plotDirectory, '/meanCNVinGroup.pdf'), genome=genome)
      plotMeanCNVtoFile(metaData, project, meanCNVs$outGroup, cosmicDirectory=cosmicDirectory, forceRedo=forceRedoMeanPlot,
                        plotFile=paste0(plotDirectory, '/meanCNVoutGroup.pdf'), genome=genome)
      
      catLog('mutation matrix...')
      plotSubgroupMutationMatrix(metaData, meanCNVs, project, subgroup, cosmicDirectory=cosmicDirectory,
                                 nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot)    
      plotMultipageMutationMatrix(metaData, meanCNVs$inGroup, project, cosmicDirectory=cosmicDirectory,
                                  nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot,
                                  plotFile=paste0(plotDirectory, '/mutationMatrixInGroup.pdf'), genome=genome)
      plotMultipageMutationMatrix(metaData, meanCNVs$outGroup, project, cosmicDirectory=cosmicDirectory,
                                  nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot,
                                  plotFile=paste0(plotDirectory, '/mutationMatrixOutGroup.pdf'), genome=genome)
      catLog('done.\n')
    }
  }
  
  catLog('Done! Invisible returning mean CNV data.\n\n')
  invisible(meanCNV)
}

#checks the mutation rates for gain, amplification, loss, complete loss, SNVs and biallelic loss
#returns matrices for genes over all samples, as well as rates in terms of fraction of individuals
#that have the gene mutated.
getMeanCNV = function(metaData, samples, variants, genome=genome, clonalityCut=0.4, cpus=1, ignoreCNAonly=F) {
  individuals = metaData$samples[samples,]$INDIVIDUAL
  uInd = unique(individuals)
  indSamples = lapply(uInd, function(i) which(individuals == i))
  cnvs = getCNV(metaData, samples)
    
  #fits = getFit(metaData, samples)
  sex = ifelse(sapply(1:length(cnvs), function(i) mean(cnvs[[i]]$clusters[xToChr(cnvs[[i]]$clusters$x1, genome)=='X',]$M)-mean(cnvs[[i]]$clusters[xToChr(cnvs[[i]]$clusters$x1, genome)=='Y',]$M)) > 3, 'female', 'male')
  if ( any(is.na(sex)) ) sex[is.na(sex)] = 'female'
  names(sex) = names(cnvs)
    
  cnvRates = getCNVrates(metaData, cnvs, sex, individuals, clonalityCut=clonalityCut, cpus=cpus, genome=genome)
    
  varSamples = names(variants$variants)
  varInd = metaData$samples[varSamples,]$INDIVIDUAL
  uVarInd = unique(varInd)
    
  snvRates = getSNVrates(metaData, variants, clonalityCut=clonalityCut, genome=genome, cpus=cpus)
  mutationRate = snvRates$mutationRate
  hitMx = snvRates$hitMx
  
  doubleLossRates = mergeDoubleLoss(snvRates, cnvRates, individuals, ignoreCNAonly=ignoreCNAonly)
  doubleLossRate = doubleLossRates$doubleLossRate
  doubleLossMx = doubleLossRates$doubleLossMx
  
  meanCNV = list(cnvs=cnvs, cnvRates=cnvRates, snvRates=snvRates, doubleLossRates=doubleLossRates, individuals=individuals)
  return(meanCNV)
}

#plots the mutation rates over the genome (CNVs and SNVs) to file.
#also outputs the data from the analysis to a spreadsheet.
plotMeanCNVtoFile = function(metaData, project, meanCNV, cosmicDirectory='', plotFile='', forceRedo=F, genome='hg19') {
  if ( plotFile == '' )
    plotFile = paste0(metaData$project[project,]$plotDirectory, '/meanCNV.pdf')
  if ( file.exists(plotFile) & !forceRedo ) return()

  pdf(plotFile, width=20, height=10)
  output = plotMeanCNV(metaData, meanCNV=meanCNV, cosmicDirectory=cosmicDirectory, add=F, genome=genome)
  dev.off()

  outputMeanCNV(metaData, project, output, outName=gsub('.pdf$', '',plotFile))
}

#plots the mutation reates over the genome.
plotMeanCNV = function(metaData, meanCNV, cosmicDirectory='', add=F, printGeneNames=T, meanCNV2=NA, genome='hg19', filterIG=T, filterTR=T, filterHLA=T, dontCountRepeatedSNVs=F, plotSex=T, outputData=T, plotMutations=T, reset=T, addLegend=T, chr.col=mcri('green', 0.3), label.cex=1, chr.cex=1, chr.stagger=0, chr.staggerFrom=1, ...) {
  samples = names(meanCNV$cnvs)
  individuals = metaData$samples[samples,]$INDIVIDUAL
  uInd = unique(individuals)
  indSamples = lapply(uInd, function(ind) which(individuals == ind))

  cL = chrLengths(genome=genome)
  if ( !plotSex ) cL = cL[!(names(cL) %in% c('X', 'Y'))]
  CR = meanCNV$cnvs[[1]]$CR
  if ( !plotSex ) CR = CR[!(xToChr(CR$x1, genome=genome) %in% c('X', 'Y')),]
  
  maxCNV = max(meanCNV$cnvRates$gain, meanCNV$cnvRates$loss, meanCNV$cnvRates$loh)
  maxLOH = max(meanCNV$cnvRates$loh)
  if ( !plotSex ) {
    use = !(xToChr(CR$x1, genome=genome) %in% c('X', 'Y'))
    maxCNV = max(meanCNV$cnvRates$gain[use], meanCNV$cnvRates$loss[use], meanCNV$cnvRates$loh[use])
    maxLOH = max(meanCNV$cnvRates$loh[use])
  }
  if ( class(meanCNV2) != 'logical' && !is.na(meanCNV2) ) {
    maxCNV = max(maxCNV, meanCNV2$cnvRates$gain, meanCNV2$cnvRates$loss, meanCNV2$cnvRates$loh)
    maxLOH = max(maxLOH, meanCNV2$cnvRates$loh)
  }
  
  spaceForSNVs = 0.5
  if ( !plotMutations ) spaceForSNVs = maxCNV*0.02
  spaceForLoh = pmin(maxLOH, maxCNV)/maxCNV
  
  indTicks = (1:length(uInd))/length(uInd)
  indTicks = indTicks[indTicks <= maxCNV+0.05]
  heavy = rep(T, length(indTicks))
  while ( sum(heavy) > 10 ) heavy[heavy] = rep(c(F,T), length.out = sum(heavy))
  cnvScale = 1/max(indTicks)
  yTicks = indTicks*cnvScale

  indTicksLOH = (1:length(uInd))/length(uInd)
  indTicksLOH = indTicksLOH[indTicksLOH <= maxLOH+0.05]
  heavyLOH = rep(T, length(indTicksLOH))
  while ( sum(heavyLOH) > 10*spaceForLoh ) heavyLOH[heavyLOH] = rep(c(F,T), length.out = sum(heavyLOH))
  cnvScaleLOH = 1/max(indTicksLOH)
  yTicksLOH = indTicksLOH*cnvScaleLOH
  
  yTicks = c(-spaceForLoh - spaceForSNVs - yTicks, -spaceForSNVs - yTicksLOH*spaceForLoh*0.85, spaceForSNVs + yTicks)
  indTicks = c(indTicks, indTicksLOH, indTicks)
  heavy = c(heavy, heavyLOH, heavy)
  ymax = max(yTicks)
  ymin = min(yTicks)

  xs = sort(c((CR$x1+CR$x2)/2, 0, cumsum(cL)))
  names(xs) = xs
  xsOri = c((CR$x1+CR$x2)/2)
  delimiter = xs %in% c(0, cumsum(cL))
  meanCNVori = meanCNV

  #split up delimiter points into two points, first matching previous value, second matching following.
  for ( cnv in c('gain', 'loss', 'amp', 'cl', 'loh') ) {
    vec = meanCNV$cnvRates[[cnv]]
    is = which(delimiter)
    vec[is] = vec[pmax(1,(is-1))]
    newIs = is + 0.5
    newValues = vec[pmin(length(vec), is+1)]
    vec = c(vec, newValues)[order(c(seq_along(vec), newIs))]
    meanCNV$cnvRates[[cnv]] = vec
  }
  xs = sort(c((CR$x1+CR$x2)/2, rep(c(0, cumsum(cL)),2)))
  Nx = length(xs)
  
  #set up plot
  if ( !add ) {
    if ( reset ) {
      par(oma=c(0,0,0,0))
      par(mar=c(0,0,0,0))
    }
    plot(xs, xs, type='n', xlim=c(-0.02*max(xs), max(xs)*(1 +0.2*addLegend)), ylim=c(ymin, ymax+(ymax-ymin)*0.03)*0.97,
         xaxt='n', yaxt='n', frame=F, xlab='', ylab='', ...)

    #horizontal lines for number of patients
    catLog('plotting CNVs..')
    segments(0, yTicks, max(xs), yTicks, lwd=1, col=rgb(0,0,0,0.05+heavy*0.15))
    text(-max(xs)*0.02, yTicks[heavy], round(indTicks[heavy], 2), col=rgb(0,0,0,1), cex=label.cex)
    text(-max(xs)*0.05, -spaceForLoh/2, 'frequency', srt=90, cex=label.cex)

    
    #mark chromosome boundaries
    superFreq:::addChromosomeLines(c(ymin, ymax*1.05), col=chr.col, lwd=1.5, genome=genome,
                                   onlyNumbers=!plotSex, cex=chr.cex, stagger=chr.stagger, staggerFrom=chr.staggerFrom)
  }

  #polygons for gain/amplification/loss/complete loss
  if ( !add ) {
    polygon(c(xs, rev(xs)), spaceForSNVs+c(meanCNV$cnvRates$gain[1:Nx], rep(0, length(xs)))*cnvScale,
            col=mcri('orange'), border=mcri('orange'))
    polygon(c(xs, rev(xs)), spaceForSNVs+c(meanCNV$cnvRates$amp[1:Nx], rep(0, length(xs)))*cnvScale,
            col=mcri('red'), border=mcri('red'))
    polygon(c(xs, rev(xs)), -spaceForSNVs-c(meanCNV$cnvRates$loh[1:Nx], rep(0, length(xs)))*cnvScaleLOH*spaceForLoh*0.85,
            col=mcri('green'), border=mcri('green'))
    polygon(c(xs, rev(xs)), -spaceForSNVs-spaceForLoh-c(meanCNV$cnvRates$loss[1:Nx], rep(0, length(xs)))*cnvScale,
            col=mcri('cyan'), border=mcri('cyan'))
    polygon(c(xs, rev(xs)), -spaceForSNVs-spaceForLoh-c(meanCNV$cnvRates$cl[1:Nx], rep(0, length(xs)))*cnvScale,
            col=mcri('blue'), border=mcri('blue'))
  } else {
    lines(xs, spaceForSNVs + meanCNV$cnvRates$gain[1:Nx]*cnvScale, col='white', lwd=2.5)
    lines(xs, spaceForSNVs + meanCNV$cnvRates$gain[1:Nx]*cnvScale, col=mcri('red'), lwd=1)
    lines(xs, spaceForSNVs + meanCNV$cnvRates$amp[1:Nx]*cnvScale, col='white', lwd=2.5)
    lines(xs, spaceForSNVs + meanCNV$cnvRates$amp[1:Nx]*cnvScale, col=mcri('darkred'), lwd=1)
    lines(xs, -spaceForSNVs - meanCNV$cnvRates$loh[1:Nx]*cnvScale*spaceForLoh, col='white', lwd=2.5)
    lines(xs, -spaceForSNVs - meanCNV$cnvRates$loh[1:Nx]*cnvScale*spaceForLoh, col=mcri('green'), lwd=1)
    lines(xs, -spaceForSNVs - spaceForLoh - meanCNV$cnvRates$loss[1:Nx]*cnvScale, col='white', lwd=2.5)
    lines(xs, -spaceForSNVs - spaceForLoh - meanCNV$cnvRates$loss[1:Nx]*cnvScale, col=mcri('blue'), lwd=1)
    lines(xs, -spaceForSNVs - spaceForLoh - meanCNV$cnvRates$cl[1:Nx]*cnvScale, col='white', lwd=2.5)
    lines(xs, -spaceForSNVs - spaceForLoh - meanCNV$cnvRates$cl[1:Nx]*cnvScale, col=mcri('darkblue'), lwd=1)
  }

  #legends
  if ( !add & addLegend ) {
    legend('bottomright', c('gain', 'amplification', 'CNN LOH', 'loss', 'complete loss'), lwd=rep(10, 5),
           col=mcri(c('orange', 'red', 'green', 'cyan', 'blue')), bg='white')
  }

  spaceForGenes = 0.15*spaceForSNVs

  if ( plotMutations ) {
    catLog('plotting mutations..')
    maxMut = max(c(0, meanCNV$snvRates$mutationRate))
    maxDL = max(meanCNV$doubleLossRates$doubleLossRate)
    
    maxSNV = max(maxMut, maxDL)
    if ( class(meanCNV2) != 'logical' && !is.na(meanCNV2) ) {
      maxSNV = max(maxSNV, meanCNV2$snvRates$mutationRate, meanCNV2$doubleLossRates$doubleLossRate)
    }
    
    indTicks = (1:length(uInd))/length(uInd)
    indTicks = indTicks[indTicks <= maxSNV]
    heavy = rep(T, length(indTicks))
    while ( sum(heavy) > 5 ) heavy[heavy] = rep(c(F,T), length.out = sum(heavy))
    snvScale = 1/max(c(0,indTicks))*(spaceForSNVs - spaceForGenes - 0.2*spaceForSNVs)
    yTicks = spaceForGenes + indTicks*snvScale
    yTicks = c(-yTicks, yTicks)
    indTicks = c(-indTicks, indTicks)
    
    if ( !add ) {
      segments(rep(0, length(yTicks)), yTicks, rep(max(xs), length(yTicks)), yTicks, lwd=1, col=mcri('cyan', 0.05+heavy*0.15))
    if ( length(indTicks[heavy]) > 0 )
      text(-max(xs)*0.02, yTicks[heavy], round(abs(indTicks[heavy]), 2), col=mcri('cyan'))
    }
    
    meanX = (CR$x1 + CR$x2)/2
    Nx = length(meanX)
    names(meanX) = rownames(meanCNV$cnvRates$gainMx)[1:Nx][!delimiter]
    SNVgenesX = meanCNV$snvRates$hitGenesX[1:Nx]
    DLgenesX = meanX[names(meanCNV$doubleLossRates$doubleLossRate)][1:Nx]
    
    if ( !printGeneNames ) {
      if ( length(SNVgenesX) > 0 )
        segments(SNVgenesX, spaceForGenes+meanCNV$snvRates$mutationRate[1:Nx]*snvScale,
                 SNVgenesX, spaceForGenes, col=mcri('cyan'), lwd=4)
      if ( length(DLgenesX) > 0 )
        segments(DLgenesX, -spaceForGenes - meanCNV$doubleLossRates$doubleLossRate[1:Nx]*snvScale,
                 DLgenesX, -spaceForGenes, col=mcri('orange'), lwd=4)
    } else {
      if ( length(SNVgenesX) > 0 )
        segments(SNVgenesX, spaceForGenes+meanCNV$snvRates$mutationRate[1:Nx]*snvScale,
                 SNVgenesX, spaceForGenes, col=mcri('blue'), lwd=1)
      if ( length(DLgenesX) > 0 )
        segments(DLgenesX, -spaceForGenes - meanCNV$doubleLossRates$doubleLossRate[1:Nx]*snvScale,
                 DLgenesX, -spaceForGenes, col=mcri('red'), lwd=1)
    }
    
    if ( printGeneNames ) {
      mutRate = meanCNV$snvRates$mutationRate[1:Nx]
      isIG = grepl('IG.[VLJC][0-9].*', names(mutRate)) | grepl('IGH[ADG][0-9].*', names(mutRate)) | grepl('LILR.*', names(mutRate))
      isTR = grepl('TR[AD][VJ][0-9].?', names(mutRate))
      isHLA = grepl('HLA-.*', names(mutRate))
      
      if ( filterIG ) mutRate = mutRate[!isIG]
      if ( filterTR ) mutRate = mutRate[!isTR]
      if ( filterHLA ) mutRate = mutRate[!isHLA]
      
      topMutGenes = names(sort(mutRate, decreasing=T)[1:min(10, length(mutRate))])
      mark = names(meanCNV$snvRates$mutationRate)[1:Nx] %in% topMutGenes
      if ( length(mark) > 0 && any(is.na(SNVgenesX[mark])) ) {
        warning('Removing genes from top 10 SNVs, as some are not found. Likely caused by unmatching gene names.')
        mark[is.na(SNVgenesX)] = F
      }
      if ( cosmicDirectory != '' )
        COSMICgenes = names(getCosmicCensusDensity(cosmicDirectory=cosmicDirectory))
      else
        COSMICgenes = c()
      printGenes = names(meanCNV$snvRates$mutationRate)[mark]
      if ( length(printGenes) > 0 ) {
        text(spreadPositions(SNVgenesX[mark], max(xs)*0.035), spaceForGenes*0.4,
             printGenes, col=ifelse(printGenes %in% COSMICgenes, mcri('green'), mcri('darkblue')),
             font=ifelse(printGenes %in% COSMICgenes, 2, 1), cex=0.6)
        segments(spreadPositions(SNVgenesX[mark], max(xs)*0.035), spaceForGenes*0.6,
                 SNVgenesX[mark], spaceForGenes, col=mcri('darkblue'))
      }
      
      dlRate = meanCNV$doubleLossRates$doubleLossRate[1:Nx]
      isIG = grepl('IG.[VLJC][0-9].*', names(dlRate)) | grepl('IGH[ADG][0-9].*', names(dlRate))
      isTR = grepl('TR[AD][VJ][0-9].?', names(dlRate))
      isHLA = grepl('HLA-.*', names(dlRate))
      
      if ( filterIG ) dlRate = dlRate[!isIG]
      if ( filterTR ) dlRate = dlRate[!isTR]
      if ( filterHLA ) dlRate = dlRate[!isHLA]
      
      topDHGenes = names(sort(dlRate, decreasing=T)[1:min(10, length(dlRate))])
      topDHGenes = topDHGenes[topDHGenes != '?']
      topDHGenes = topDHGenes[meanCNV$doubleLossRates$doubleLossRate[1:Nx][topDHGenes] > 0]
      mark = names(meanCNV$doubleLossRates$doubleLossRate)[1:Nx] %in% topDHGenes
      if ( any(is.na(DLgenesX[mark])) ) {
        warning('Removing genes from top 10 SNVs, as some are not found. Likely caused by unmatching gene names.')
        mark[is.na(DLgenesX)] = F
      }
      printGenes = names(meanCNV$doubleLossRates$doubleLossRate)[1:Nx][mark]
      if ( length(printGenes) > 0 ) {
        text(spreadPositions(DLgenesX[mark], max(xs)*0.035), -spaceForGenes*0.4,
             printGenes, col=ifelse(printGenes %in% COSMICgenes, mcri('green'), mcri('red')),
             font=ifelse(printGenes %in% COSMICgenes, 2, 1), cex=0.6)
        segments(spreadPositions(DLgenesX[mark], max(xs)*0.035), -spaceForGenes*0.6,
                 DLgenesX[mark], -spaceForGenes, col=mcri('red'))
      }
    }
    
    if ( !add ) {
      legend('right', c('changing SNV', 'biallelic loss'), col=mcri(c('darkblue', 'red')), lwd=3)
    }
    catLog('done.\n')
  }

  if ( outputData ) {
    catLog('Setting up output...')
    meanCNV = meanCNVori
    xs = sort(c((meanCNV$cnvs[[1]]$CR$x1+meanCNV$cnvs[[1]]$CR$x2)/2, 0, cumsum(chrLengths(genome=genome))))
    genes = rownames(meanCNV$cnvRates$gainMx)[!delimiter]
    allMutationRates = rep(0, length(genes))
    names(allMutationRates) = genes
    use = names(meanCNV$snvRates$mutationRate) %in% names(allMutationRates)
    allMutationRates[names(meanCNV$snvRates$mutationRate)[use]] = meanCNV$snvRates$mutationRate[use]
    allDLRates = rep(0, length(genes))
    names(allDLRates) = genes
    use = names(meanCNV$doubleLossRates$doubleLossRate) %in% names(allDLRates)
    allDLRates[names(meanCNV$doubleLossRates$doubleLossRate)[use]] = meanCNV$doubleLossRates$doubleLossRate[use]
    
    overview = data.frame(chr=xToChr(xs[!delimiter], genome), pos0=round(xToPos(xs[!delimiter], genome)), gene=genes,
      gainRate=meanCNV$cnvRates$gain[!delimiter], ampRate=meanCNV$cnvRates$amp[!delimiter],
      lossRate=meanCNV$cnvRates$loss[!delimiter], clRate=meanCNV$cnvRates$cl[!delimiter],
      lohRate=meanCNV$cnvRates$loh[!delimiter], mutationRate=allMutationRates, doubleLossRate=allDLRates)
    
    gains = data.frame(chr=xToChr(xs[!delimiter], genome), pos=xToPos(xs[!delimiter], genome), gene=genes, meanCNV$cnvRates$gainMx[!delimiter,])
    losss = data.frame(chr=xToChr(xs[!delimiter], genome), pos=xToPos(xs[!delimiter], genome), gene=genes, meanCNV$cnvRates$lossMx[!delimiter,])
    amps = data.frame(chr=xToChr(xs[!delimiter], genome), pos=xToPos(xs[!delimiter], genome), gene=genes, meanCNV$cnvRates$ampMx[!delimiter,])
    cls = data.frame(chr=xToChr(xs[!delimiter], genome), pos=xToPos(xs[!delimiter], genome), gene=genes, meanCNV$cnvRates$clMx[!delimiter,])
    lohs = data.frame(chr=xToChr(xs[!delimiter], genome), pos=xToPos(xs[!delimiter], genome), gene=genes, meanCNV$cnvRates$lohMx[!delimiter,])
  mutations = data.frame(gene=rownames(meanCNV$snvRates$hitMx), meanCNV$snvRates$hitMx)
    DLuse = meanCNV$doubleLossRates$doubleLossRate > 0 & !duplicated(names(meanCNV$doubleLossRates$doubleLossRate))
    DLs = data.frame(gene=names(meanCNV$doubleLossRates$doubleLossRate[DLuse]),
      meanCNV$doubleLossRates$doubleLossMx[DLuse,])
    
    output = list(overview=overview, gain=gains, loss=losss, amp=amps, lohs=lohs, completeLoss=cls, mutations=mutations, doubleLoss=DLs)
    catLog('done.\n')
    invisible(output)
  }
}


outputMeanCNV = function(metaData, project, output, subgroup='', inGroup=T, outName='') {
  plotDirectory = metaData$project[project,]$plotDirectory
  if ( subgroup != '' ) plotDirectory = paste0(metaData$project[project,]$plotDirectory, '/', subgroup)
  ensureDirectoryExists(plotDirectory)

  catLog('printing to csv...')
  csvFile = paste0(plotDirectory, '/meanCNV.csv')
  if ( !inGroup ) csvFile = paste0(plotDirectory, '/meanCNVnotInGroup.csv')
  if ( outName != '' ) csvFile = paste0(outName, '.csv')
  write.csv(output$overview, csvFile)
  output = lapply(output, function(dat) dat = dat[1:min(65000, nrow(dat)), 1:min(255, ncol(dat))])
  catLog('printing to xls...')
  xlsFile = paste0(plotDirectory, '/meanCNV.xls')
  if ( !inGroup ) xlsFile = paste0(plotDirectory, '/meanCNVnotInGroup.xls')
  if ( outName != '' ) xlsFile = paste0(outName, '.xls')
  WriteXLS('output', xlsFile)
    
  catLog('done.\n')
}


getCNVrates = function(metaData, cnvs, sex, individuals, clonalityCut=0.4, cpus=1, genome='hg19') {
  uInd = unique(individuals)
  indSamples = lapply(uInd, function(i) which(individuals == i))
  
  xs = sort(c((cnvs[[1]]$CR$x1+cnvs[[1]]$CR$x2)/2, 0, cumsum(chrLengths(genome=genome))))
  names(xs) = xs
  chr = xToChr(xs, genome)

  catLog('Calculating gain..')
  gainList = mclapply(cnvs, function(cnv) sapply(xs, function(x) {
    i = which(cnv$clusters$x1 < x & cnv$clusters$x2 > x)[1]
    return(((cnv$cluster[i,]$call %in% c('AAB', 'AAA') & cnv$cluster[i,]$clonality > clonalityCut) | cnv$cluster[i,]$M > log2(1 + clonalityCut/2))[1])
    }), mc.cores=cpus)
  gainList = lapply(gainList, function(gain) ifelse(is.na(gain), F, gain))
  gainMx = do.call(cbind, gainList)
  indGain = apply(gainMx, 1, FUN=function(hasGain) sum(sapply(indSamples, function(is) any(hasGain[is]))))
  gain = indGain/length(uInd)

  catLog('loss..')
  loss = mclapply(names(cnvs), function(name) sapply(xs, function(x) {
    cnv = cnvs[[name]]
    i = which(cnv$clusters$x1 < x & cnv$clusters$x2 > x)[1]
    return(((cnv$cluster[i,]$call %in% c('A', 'CL') & cnv$cluster[i,]$clonality > clonalityCut) | cnv$cluster[i,]$M < log2(1 - clonalityCut/2))[1] & !((xToChr(x, genome) %in% c('X','Y') & sex[name] == 'male') | xToChr(x, genome) == 'Y'))
  }), mc.cores=cpus)
  lossList = lapply(loss, function(loss) ifelse(is.na(loss), F, loss))
  lossMx = do.call(cbind, lossList)
  indLoss = apply(lossMx, 1, FUN=function(hasLoss) sum(sapply(indSamples, function(is) any(hasLoss[is]))))
  loss = indLoss/length(uInd)
  
  catLog('amplification..')
  amp = mclapply(cnvs, function(cnv) sapply(xs, function(x) {
    i = which(cnv$clusters$x1 < x & cnv$clusters$x2 > x)[1]
    return(((cnv$cluster[i,]$call %in% c('AAAB', 'AAAAB', 'AAAAAB', 'AAAAAAB') & cnv$cluster[i,]$clonality > clonalityCut) | cnv$cluster[i,]$M > log2(1.6))[1])
    }), mc.cores=cpus)
  ampList = lapply(amp, function(amp) ifelse(is.na(amp), F, amp))
  ampMx = do.call(cbind, ampList)
  indAmp = apply(ampMx, 1, FUN=function(hasAmp) sum(sapply(indSamples, function(is) any(hasAmp[is]))))
  amp = indAmp/length(uInd)
  
  catLog('complete loss..')
  cl = mclapply(names(cnvs), function(name) sapply(xs, function(x) {
    cnv = cnvs[[name]]
    i = which(cnv$clusters$x1 < x & cnv$clusters$x2 > x)[1]
    return(((cnv$cluster[i,]$call == 'CL' & cnv$cluster[i,]$clonality > clonalityCut) | cnv$cluster[i,]$M < log2(0.4))[1] & !((xToChr(x, genome) == 'X' & sex[name] == 'male') | xToChr(x, genome) == 'Y'))
  }), mc.cores=cpus)
  clList = lapply(cl, function(cl) ifelse(is.na(cl), F, cl))
  clMx = do.call(cbind, clList)
  indCl = apply(clMx, 1, FUN=function(hasCl) sum(sapply(indSamples, function(is) any(hasCl[is]))))
  cl = indCl/length(uInd)

  catLog('CNN LOH..')
  loh = mclapply(names(cnvs), function(name) sapply(xs, function(x) {
    cnv = cnvs[[name]]
    i = which(cnv$clusters$x1 < x & cnv$clusters$x2 > x)[1]
    return(((cnv$cluster[i,]$call == 'AA' & cnv$cluster[i,]$clonality > clonalityCut) | (abs(cnv$cluster[i,]$M) < 0.1 & cnv$cluster[i,]$f < 0.5 - clonalityCut/2))[1] & !((xToChr(x, genome) == 'X' & sex[name] == 'male') | xToChr(x, genome) == 'Y'))
  }), mc.cores=cpus)
  lohList = lapply(loh, function(loh) ifelse(is.na(loh), F, loh))
  lohMx = do.call(cbind, lohList)
  indLoh = apply(lohMx, 1, FUN=function(hasLoh) sum(sapply(indSamples, function(is) any(hasLoh[is]))))
  loh = indLoh/length(uInd)

  
  genes = rownames(cnvs[[1]]$CR)
  delimiters = c(0, cumsum(chrLengths(genome=genome)))
  rownames = c(genes, delimiters)[order(c(cnvs[[1]]$CR$x1, delimiters))]
  colnames(gainMx) = colnames(ampMx) = colnames(lossMx) = colnames(clMx) = colnames(lohMx) = names(cnvs)
  rownames(gainMx) = rownames(ampMx) = rownames(lossMx) = rownames(clMx) = rownames(lohMx) = rownames

  return(list(gain=gain, amp=amp, loss=loss, cl=cl, loh=loh, gainMx=gainMx, ampMx=ampMx, lossMx=lossMx, clMx=clMx, lohMx=lohMx))
}


getSNVrates = function(metaData, variants, clonalityCut=0.4, genome='hg19', cpus=1) {
  variants = cleanVariants(variants)
  SNPs = variants$SNPs
  samples = names(variants$variants)
  individuals = metaData$samples[samples,]$INDIVIDUAL
  uInd = unique(individuals)
  indSamples = lapply(uInd, function(i) which(individuals == i))

  catLog('somatic SNVs')
  hitX = mclapply(variants$variants, function(q) {
    catLog('.')
    q = q[q$somaticP > 0.1 & q$var > q$cov*clonalityCut/2 & q$severity <= 11,]
    return(list('x'=q$x, 'rowname'=paste0(q$x, q$variant), 'sev'=q$severity, 'genes'=q$inGene))
  }, mc.cores=cpus)
  hitInfo = data.frame(x=unlist(lapply(hitX, function(hit) hit$x)), rowname=unlist(lapply(hitX, function(hit) hit$rowname)), severity=unlist(lapply(hitX, function(hit) hit$sev)), genes=unlist(lapply(hitX, function(hit) hit$genes)), stringsAsFactors=F)
  hitInfo = hitInfo[!duplicated(hitInfo$rowname),]
  xMx = do.call(cbind, lapply(hitX, function(hit) hitInfo$rowname %in% hit$rowname))
  rownames(xMx) =rownames(hitInfo) = hitInfo$rowname
  xRate = rowsums(xMx)
  
  catLog('mutated genes..')
  hitGenes = mclapply(variants$variants, function(q) {
    catLog('.')
    q = q[q$somaticP > 0.1 & q$var > q$cov*clonalityCut/2 & q$severity <= 11,]
    if ( nrow(q) == 0 ) return(c())
    return(q$inGene)
  }, mc.cores=cpus)
  catLog('\n')
  mutatedGenes = unique(unlist(lapply(hitGenes, unique)))
  mutatedGenes = mutatedGenes[mutatedGenes != '?']
  hitMx = do.call(cbind, lapply(hitGenes, function(genes) mutatedGenes %in% genes))
  rownames(hitMx) = mutatedGenes
  indHit = apply(hitMx, 1, FUN=function(hasHit) sum(sapply(indSamples, function(is) any(hasHit[is]))))
  mutationRate = indHit/length(uInd)

  doubleHitGenes = unique(unlist(lapply(hitGenes, function(genes) genes[duplicated(genes)])))
  doubleHitGenes = doubleHitGenes[doubleHitGenes != '?']
  doubleHitMx = do.call(cbind, lapply(hitGenes, function(genes) doubleHitGenes %in% genes[duplicated(genes)]))
  rownames(doubleHitMx) = doubleHitGenes
  indDoubleHit = apply(doubleHitMx, 1, FUN=function(hasHit) sum(sapply(indSamples, function(is) any(hasHit[is]))))
  doubleMutationRate = indDoubleHit/length(uInd)

  hitSNPs = SNPs[SNPs$inGene %in% mutatedGenes,]
  hitGenesX = unlist(mclapply(mutatedGenes, function(gene) mean(hitSNPs[hitSNPs$inGene == gene,]$x), mc.cores=cpus))

  return(list(hitGenesX=hitGenesX, mutationRate=mutationRate, hitMx=hitMx, xRate=xRate, xMx=xMx, hitInfo=hitInfo,
              doubleHitMx=doubleHitMx, doubleMutationRate=doubleMutationRate))
}


mergeDoubleLoss = function(snvRates, cnvRates, individuals, ignoreCNAonly=F) {
  lossMx = cnvRates$lossMx
  clMx = cnvRates$clMx
  lohMx = cnvRates$lohMx
  hitMx = snvRates$hitMx
  doubleHitMx = snvRates$doubleHitMx

  if ( !setequal(colnames(lossMx), colnames(hitMx)) ) {
    warning('Couldnt do merged snv hits for the double loss analysis, as snv and cnv samples didnt match.')
    print(colnames(lossMx))
    print(colnames(hitMx))
    return(list('doubleLossMx'=clMx, 'doubleLossRate'=cnvRates$cl))
  }

  if ( !ignoreCNAonly ) lossMx = lossMx + clMx
  hitUse = rownames(hitMx) %in% rownames(lossMx)
  lossMx[rownames(hitMx[hitUse,]),] = lossMx[rownames(hitMx[hitUse,]),] + hitMx[hitUse,colnames(lossMx)]
  lossMx[rownames(hitMx[hitUse,]),] = lossMx[rownames(hitMx[hitUse,]),] + hitMx[hitUse,colnames(lossMx)]*lohMx[rownames(hitMx[hitUse,]),]
  doubleHitUse = rownames(doubleHitMx) %in% rownames(lossMx)
  lossMx[rownames(doubleHitMx[doubleHitUse,]),] = lossMx[rownames(doubleHitMx[doubleHitUse,]),] + doubleHitMx[doubleHitUse,colnames(doubleHitMx)]

  doubleLossMx = lossMx > 1
  rownames(doubleLossMx) = rownames(lossMx)

  uInd = unique(individuals)
  indSamples = lapply(uInd, function(i) which(individuals == i))
  indDoubleLoss = apply(doubleLossMx, 1, FUN=function(hasDL) sum(sapply(indSamples, function(is) any(hasDL[is]))))
  doubleLossRate = indDoubleLoss/length(uInd)
  names(doubleLossRate) = rownames(lossMx)

  return(list('doubleLossMx'=doubleLossMx, 'doubleLossRate'=doubleLossRate))
}


plotMultipageMutationMatrix =
  function(metaData, meanCNV, project, meanCNVs=NA, cosmicDirectory='',
           priorCnvWeight=1, nGenes=30, dontCountRepeatedSNVs=F,
           pages=1, forceRedo=F, plotFile='', cancerGenes=NA, genome='hg19') {
  if ( plotFile == '' )
    plotFile = paste0(metaData$project[project,]$plotDirectory, '/mutationMatrix.pdf')
  if ( file.exists(plotFile) & !forceRedo ) return()

  pdf(plotFile, width=20, height=10)
  for ( page in 1:pages ) {
    plotMutationMatrix(metaData, meanCNV, cosmicDirectory=cosmicDirectory, priorCnvWeight=priorCnvWeight, dontCountRepeatedSNVs=dontCountRepeatedSNVs,
                       nGenes=nGenes, skipFirst=(page-1)*nGenes, forceRedo=forceRedo, cancerGenes=cancerGenes, genome=genome)
  }
  dev.off()
  
}

cohortGeneScore = function(metaData, meanCNV, priorCnvWeight=1, filterIG=T, filterTR=T, filterHLA=T, genome='hg19', dontCountRepeatedSNVs=F) {
    
  xs = sort(c((meanCNV$cnvs[[1]]$CR$x1+meanCNV$cnvs[[1]]$CR$x2)/2, 0, cumsum(chrLengths(genome=genome))))
  delimiter = xs %in% c(0, cumsum(chrLengths(genome=genome)))
  
  samples = colnames(meanCNV$cnvRates$gainMx)
  individuals = sampleToIndividual(metaData, metaData$samples$NAME)
  normals = metaData$samples$NORMAL
  names(individuals) = names(normals) = metaData$samples$NAME
  hasMatchedNormal = !is.na(findCorrespondingNormal(metaData$samples$NAME, individuals, normals))[samples]
    
  weightedxMx = t((1+hasMatchedNormal)*t(meanCNV$snvRates$xMx))
  if ( dontCountRepeatedSNVs ) {
    weightedxMx = t(apply(weightedxMx, 1, function(row) {
      row[1:length(row) != which(row == max(row))[1]] = 0
      return(row)
      }))
  }
  individuals = sampleToIndividual(metaData, samples)
  uInd = unique(individuals)
  indSamples = lapply(uInd, function(i) which(individuals == i))
  snvGenes = unique(meanCNV$snvRates$hitInfo$genes)
  if ( length(snvGenes) == 0 ) return(c())
  geneMx = lapply(snvGenes, function(gene) weightedxMx[meanCNV$snvRates$hitInfo$genes==gene,,drop=F])
  indMut = sapply(geneMx, function(mx) sum(sapply(indSamples, function(is) {
    indMx = mx[,is,drop=F]
    return(min(2,sum(rowSums(indMx) > 0)))
  })))
  weightedMutationRate = indMut/length(uInd)
  names(weightedMutationRate) = snvGenes


  #weightedHitMx = t((1+hasMatchedNormal)*t(meanCNV$snvRates$hitMx))
  #individuals = sampleToIndividual(metaData, samples)
  #uInd = unique(individuals)
  #indSamples = lapply(uInd, function(i) which(individuals == i))
  #indMut = apply(weightedHitMx, 1, FUN=function(weight) sum(sapply(indSamples, function(is) max(weight[is]))))
  #weightedMutationRate = indMut/length(uInd)

  if ( dontCountRepeatedSNVs ) {
    names(geneMx) = snvGenes
    for ( gene in snvGenes ) {
      if ( !(gene %in% rownames(meanCNV$cnvRates$lossMx)) ) next
      lostOne = meanCNV$cnvRates$lossMx[gene,]
      notSNV = colSums(geneMx[[gene]]) == 0
      dl = meanCNV$doubleLossRates$doubleLossMx[gene,]
      fp = lostOne & notSNV & dl
      if ( any(fp) ) {
        meanCNV$doubleLossRates$doubleLossMx[gene,fp] = FALSE
        meanCNV$doubleLossRates$doubleLossRate[gene] =
          sum(meanCNV$doubleLossRates$doubleLossMx[gene,])/ncol(meanCNV$doubleLossRates$doubleLossMx)
      }
    }
  }

  genes = rownames(meanCNV$cnvRates$gainMx)[!delimiter]
  allMutationRates = rep(0, length(genes))
  names(allMutationRates) = genes
  use = names(weightedMutationRate) %in% names(allMutationRates)
  allMutationRates[names(weightedMutationRate)[use]] = weightedMutationRate[use]
  allDLRates = rep(0, length(genes))
  names(allDLRates) = genes
  use = names(meanCNV$doubleLossRates$doubleLossRate) %in% names(allDLRates)
  allDLRates[names(meanCNV$doubleLossRates$doubleLossRate)[use]] = meanCNV$doubleLossRates$doubleLossRate[use]

  gain = meanCNV$cnvRates$gain[!delimiter]
  amp = meanCNV$cnvRates$amp[!delimiter]
  loss = meanCNV$cnvRates$loss[!delimiter]
  cl = meanCNV$cnvRates$cl[!delimiter]
  loh = meanCNV$cnvRates$loh[!delimiter]

  mutation = allMutationRates
  
  doubleLoss = allDLRates
  names(gain) = names(amp) = names(loss) = names(cl) = names(loh) = names(mutation) = names(doubleLoss) = genes

  isIG = grepl('IG.[VLJC][0-9].*', genes) | grepl('IGH[ADG][0-9].*', genes) | grepl('LILR.*', genes)
  isTR = grepl('TR[AD][VJ][0-9].?', genes)
  isHLA = grepl('HLA-.*', genes)

  cnvScore = gain*0.5 + loss*0.5 + amp + cl + loh*0.5
  cnvWeight = sqrt(mean(mutation)/mean(cnvScore))*priorCnvWeight
  score = cnvScore*cnvWeight + mutation + (doubleLoss - cl + cl*cnvWeight)*0.5 #DL downweighted, as it is already counted as loss/cl/mutation
  if ( filterIG ) score = score*(1-isIG)
  if ( filterTR ) score = score*(1-isTR)
  if ( filterHLA ) score = score*(1-isHLA)

  return(score)
}


#This function takes the output from meanCNV, selects the most interesting genes
#and plots the mutation and CNV status of those genes for all samples.
plotMutationMatrix = function(metaData, meanCNV, cosmicDirectory='', priorCnvWeight=1, nGenes=30, filterIG=T, dontCountRepeatedSNVs=F,
                              skipFirst=0, forceRedo=F, meanCNV2=NA, add=F, cancerGenes=NA, genome='hg19') {
  if ( is.na(cancerGenes[1]) ) {
    score = cohortGeneScore(metaData, meanCNV, priorCnvWeight=priorCnvWeight, filterIG=filterIG, dontCountRepeatedSNVs=dontCountRepeatedSNVs)
    score[names(score)=='?'] = 0
    if ( class(meanCNV2) != 'logical' && !is.na(meanCNV2) ) {
      score2 = cohortGeneScore(metaData, meanCNV2, priorCnvWeight=priorCnvWeight, filterIG=filterIG, genome=genome, dontCountRepeatedSNVs=dontCountRepeatedSNVs)
      score = abs(score-score2)
    }
    cancerGenes = names(sort(score, decreasing=T)[skipFirst + 1:nGenes])
  }
  else {
    cancerGenes = cancerGenes[cancerGenes %in% rownames(meanCNV$cnvRates$gainMx)][skipFirst + 1:nGenes]
  }
  cancerGenes = cancerGenes[!is.na(cancerGenes)]
  if ( length(cancerGenes) == 0 ) return()

  xs = sort(c((meanCNV$cnvs[[1]]$CR$x1+meanCNV$cnvs[[1]]$CR$x2)/2, 0, cumsum(chrLengths(genome=genome))))
  delimiter = xs %in% c(0, cumsum(chrLengths(genome=genome)))

  gains = meanCNV$cnvRates$gainMx[!delimiter,,drop=F]
  losss = meanCNV$cnvRates$lossMx[!delimiter,,drop=F]
  amps = meanCNV$cnvRates$ampMx[!delimiter,,drop=F]
  cls = meanCNV$cnvRates$clMx[!delimiter,,drop=F]
  lohs = meanCNV$cnvRates$lohMx[!delimiter,,drop=F]
  mutations = gains*0
  doubleSNV = gains*0
  use = rownames(meanCNV$snvRates$hitMx) %in% rownames(gains)
  mutations[rownames(meanCNV$snvRates$hitMx)[use],] = meanCNV$snvRates$hitMx[use,,drop=F]
  use = rownames(meanCNV$snvRates$doubleHitMx) %in% rownames(gains)
  doubleSNV[rownames(meanCNV$snvRates$doubleHitMx)[use],] = meanCNV$snvRates$doubleHitMx[use,,drop=F]
  DLs = meanCNV$doubleLossRates$doubleLossMx

  samples = colnames(gains)
  samples = samples[order(metaData$samples[samples,]$INDIVIDUAL)]
  individuals = metaData$samples[samples,]$INDIVIDUAL
  gains = gains[cancerGenes,samples,drop=F]
  losss = losss[cancerGenes,samples,drop=F]
  amps = amps[cancerGenes,samples,drop=F]
  cls = cls[cancerGenes,samples,drop=F]
  lohs = lohs[cancerGenes,samples,drop=F]
  mutations = mutations[cancerGenes,samples,drop=F]
  doubleSNV = doubleSNV[cancerGenes,samples,drop=F]
  DLs = DLs[cancerGenes,samples,drop=F]

  bgCol = rgb(0.9, 0.9, 0.9)
  cnvsToCols = function(gains, amps, losss, cls, lohs) {
    is = 1:(nrow(gains)*ncol(gains))
    cols = sapply(is, function(i) cnvToCol(gains[i], amps[i], losss[i], cls[i], lohs[i]))
    return(matrix(cols, nrow=nrow(gains)))
  }
  cnvToCol = function(gain, amp, loss, cl, loh) {
    if ( cl ) return(mcri('blue'))
    if ( loss ) return(mcri('cyan'))
    if ( loh ) return(mcri('green'))
    if ( amp ) return(mcri('darkred'))
    if ( gain ) return(mcri('red'))
    return(bgCol)
  }

  cnvSize = 1
    
  x = 1:ncol(gains) + sapply(individuals, function(ind) which(unique(individuals) == ind))/3
  names(x) = samples
  y = nrow(gains):1
  names(y) = cancerGenes

  
  xmax = max(x+0.5)
  ymax = max(y+0.5)

  if ( class(meanCNV2) != 'logical' && !is.na(meanCNV2) ) {
    individuals2 = sampleToIndividual(metaData, names(meanCNV2$cnvs))
    maxX2 = ncol(meanCNV2$cnvRates$gainMx) + length(unique(individuals2))/3
    if ( add ) x = x + maxX2 + 1
    else xmax = xmax + maxX2 + 1
  }
  
  squareX = sapply(x, function(X) rep(X, length(y)))
  squareY = t(sapply(y, function(Y) rep(Y, length(x))))
  indX = sapply(unique(individuals), function(ind) mean(x[ind == individuals]))
  cnvCol = cnvsToCols(gains, amps, losss, cls, lohs)
  frameCol = matrix(ifelse(DLs, mcri('red'), rgb(0,0,0,0)), nrow=nrow(gains))
  snvCol = matrix(ifelse(mutations, mcri('red'), bgCol), nrow=nrow(gains))
  doubleSnvCol = matrix(ifelse(doubleSNV, mcri('white'), snvCol), nrow=nrow(gains))

  if ( !add ) {
    resetMargins()
    par(mar=c(0,0,0,0))
    plot(0, type='n', xlim=c(-xmax*0.02, xmax*1.2), ylim=c(-0.05*ymax,ymax*1.05), frame.plot=F, xaxt='n', yaxt='n')
  }
  rect(squareX-0.5, squareY-0.5, squareX+0.5, squareY+0.5, col=frameCol, border=F)
  rect(squareX-0.45, squareY-0.45, squareX+0.45, squareY+0.45, border=F, col=bgCol)
  if ( add ) segments(min(x)-1.25, 0.5, min(x)-1.25, ymax, lwd=5, col='black')

  boxScale = pmin(1, sqrt(20/length(indX)))
  indScale = pmin(1, sqrt(20/length(indX)))

  rect(squareX+0.25, squareY-0.25*cnvSize, squareX+0.25+0.15*cnvSize, squareY+0.25*cnvSize, col=cnvCol, border=F)
  points(squareX-0.15, squareY+0.2, pch=16, col=snvCol, cex=1.5*boxScale)
  points(squareX-0.15, squareY+0.2, pch=16, col=doubleSnvCol, cex=0.9*boxScale)

  if ( cosmicDirectory != '' ) {
    COSMICgenes = names(getCosmicCensusDensity(cosmicDirectory=cosmicDirectory))
    if ( file.exists(paste0(cosmicDirectory, '/CCGD_export.csv')) ) {
      CCGDdata = read.table(paste0(cosmicDirectory, '/CCGD_export.csv'), sep=',', header=T, stringsAsFactors=F)
      censusGenes = unique(CCGDdata$Mouse.Symbol)
      COSMICgenes = c(COSMICgenes, censusGenes)
    }
  }
  else
    COSMICgenes = c()
  if ( !add ) text(0.3, y, cancerGenes, adj=c(1, 0.5), srt=0, cex=1,
                   col=ifelse(cancerGenes %in% COSMICgenes, mcri('green'), 'black'),
                   font=ifelse(cancerGenes %in% COSMICgenes, 2, 1))
  text(indX, ymax, unique(individuals), srt=30, adj=c(0, 0), cex=indScale)
  text(x, min(y)-0.5, metaData$samples[samples,]$TIMEPOINT, srt=-30, adj=c(0, 1), cex=boxScale)

  legend('bottomright', c('SNV', '2 SNVs', 'complete loss', 'loss', 'loh', 'gain', 'amp', 'biallelic loss'), col=mcri(c('red', 'red', 'blue', 'cyan', 'green', 'red', 'darkred', 'red')), pch=c(16, 1, 15, 15, 15, 15, 15, 22), pt.cex=c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2), pt.lwd=c(1,5,1,1,1,1,3), pt.bg=mcri(c('red', 'white', 'blue', 'cyan', 'green', 'red', 'darkred', 'white')), bg='white')
 
}


#takes a meanCNV object of a project and splits it up into two meanCNV objects (returned in a list)
#one containing the samples in the provided subgroup. the other containing the samples not in the subgroup
splitMeanCNV = function(metaData, meanCNV, project, subgroup, includeNormal=F) {
  samples = names(meanCNV$cnvs)
  groupSamples = inSubgroup(metaData, project, subgroup, includeNormal=includeNormal, onlyDNA=T)
  groupSamples = samples[samples %in% groupSamples]
  otherSamples = samples[!(samples %in% groupSamples)]
  meanCNVinSubgroup = meanCNVoutsideSubgroup = meanCNV
  
  individualsIn = sampleToIndividual(metaData, groupSamples)
  meanCNVinSubgroup$individuals = individualsIn
  uIndIn = unique(individualsIn)
  indSamplesIn = lapply(uIndIn, function(i) which(individualsIn == i))
  individualsOut = sampleToIndividual(metaData, otherSamples)
  uIndOut = unique(individualsOut)
  indSamplesOut = lapply(uIndOut, function(i) which(individualsOut == i))
  meanCNVoutsideSubgroup$individuals = individualsOut

  meanCNVinSubgroup$cnvs = meanCNVinSubgroup$cnvs[groupSamples]
  meanCNVoutsideSubgroup$cnvs = meanCNVoutsideSubgroup$cnvs[otherSamples]
  entries =
    list(list(data='cnvRates', rate='gain', matrix='gainMx'),
         list(data='cnvRates', rate='amp', matrix='ampMx'),
         list(data='cnvRates', rate='loss', matrix='lossMx'),
         list(data='cnvRates', rate='cl', matrix='clMx'),
         list(data='cnvRates', rate='loh', matrix='lohMx'),
         list(data='snvRates', rate='mutationRate', matrix='hitMx'),
         list(data='snvRates', rate='xRate', matrix='xMx'),
         list(data='snvRates', rate='doubleMutationRate', matrix='doubleHitMx'),
         list(data='doubleLossRates', rate='doubleLossRate', matrix='doubleLossMx'))
  for ( entry in entries ) {
    dat = entry$data
    mx = entry$matrix
    meanCNVinSubgroup[[dat]][[mx]] = meanCNVinSubgroup[[dat]][[mx]][,groupSamples,drop=F]
    meanCNVoutsideSubgroup[[dat]][[mx]] = meanCNVoutsideSubgroup[[dat]][[mx]][,otherSamples,drop=F]
    
    individualsHitIn = apply(meanCNVinSubgroup[[dat]][[mx]], 1,
                             function(hasHit) sum(sapply(indSamplesIn, function(is) any(hasHit[is]))))
    individualsHitOut = apply(meanCNVoutsideSubgroup[[dat]][[mx]], 1,
                             function(hasHit) sum(sapply(indSamplesOut, function(is) any(hasHit[is]))))
    meanCNVinSubgroup[[dat]][[entry$rate]] = individualsHitIn/length(uIndIn)
    meanCNVoutsideSubgroup[[dat]][[entry$rate]] = individualsHitOut/length(uIndOut)
  }

  return(list('inGroup'=meanCNVinSubgroup, 'outGroup'=meanCNVoutsideSubgroup))
}


#plots the mutation rates over the genome (CNVs and SNVs) to file.
#also outputs the data from the analysis to a spreadsheet.
#this version divides the samples of the project by the provided subgroup
#and highlights differences.
plotSubgroupCNVtoFile = function(metaData, meanCNV, project, subgroup, includeNormal=F, cosmicDirectory='', add=F, forceRedo=F, genome='hg19') {
  plotDirectory = paste0(metaData$project[project,]$plotDirectory, '/', subgroup)
  ensureDirectoryExists(plotDirectory)
  plotFile = paste0(plotDirectory, '/meanCNVs.pdf')

  catLog('Splitting stats by subgroup...')
  meanCNVs = splitMeanCNV(metaData, meanCNV, project, subgroup, includeNormal)
  catLog('done.\n')

  if ( file.exists(plotFile) & !forceRedo ) return(meanCNVs)

  pdf(plotFile, width=20, height=10)
  outputOut = plotMeanCNV(metaData, meanCNV=meanCNVs$outGroup, cosmicDirectory=cosmicDirectory,
                          add=F, printGeneNames=F, meanCNV2=meanCNVs$inGroup, genome=genome)
  outputIn = plotMeanCNV(metaData, meanCNV=meanCNVs$inGroup, cosmicDirectory=cosmicDirectory,
                         add=T, printGeneNames=T, genome=genome)
  legend('topright', c(subgroup, subgroup, subgroup, paste0('not ', subgroup), paste0('not ', subgroup), paste0('not ', subgroup)), col=mcri(c('blue', 'red', 'green', 'cyan', 'orange', 'green')), lwd=c(2, 2, 2, 10, 10, 10), bg='white')
  dev.off()

  outputMeanCNV(metaData, project, outputIn, subgroup=subgroup, inGroup=T)
  outputMeanCNV(metaData, project, outputOut, subgroup=subgroup, inGroup=F)
  return(meanCNVs)
}


plotSubgroupMutationMatrix = function(metaData, meanCNVs, project, subgroup, cosmicDirectory='', priorCnvWeight=1,
                                      nGenes=30, pages=10, forceRedo=F, genome='hg19') {
  plotDirectory = paste0(metaData$project[project,]$plotDirectory, '/', subgroup)
  ensureDirectoryExists(plotDirectory)
  plotFile = paste0(plotDirectory, '/mutationMatrix.pdf')
  if ( file.exists(plotFile) & !forceRedo ) return()

  pdf(plotFile, width=20, height=10)
  for ( page in 1:pages ) {
    plotMutationMatrix(metaData, meanCNVs$inGroup, cosmicDirectory=cosmicDirectory, priorCnvWeight=priorCnvWeight,
                       nGenes=nGenes, skipFirst=(page-1)*nGenes, meanCNV2=meanCNVs$outGroup, genome=genome)
    plotMutationMatrix(metaData, meanCNVs$outGroup, cosmicDirectory=cosmicDirectory, priorCnvWeight=priorCnvWeight,
                       nGenes=nGenes, skipFirst=(page-1)*nGenes, meanCNV2=meanCNVs$inGroup, add=T, genome=genome)
  }
  dev.off()
  catLog('done.\n')
}


#' Compares sets of individuals for reccuring mutations
#'
#' @param metaDataFile character: Path to the metaData.
#' @param outputDirectories A named list of output directories, containing the entry Rdirectory and plotDirectory where the saved data and plots will be stored respectively.
#' @param project character: The project containing the subgroups.
#' @param subgroup1 character: The first subgroup(s).
#' @param subgroup2 character: The second subgroup(s).
#' @param name character: The name of the comparison. This names the output directory.
#' @param clonalityCut numeric: the minimum required clonality to be included in the analysis. Deafult 0.4.
#' @param cosmicDirectory character: The directory with the COSMIC data.
#' @param cpus integer: The maximum number of parallel processes.
#' @param forceRedoMean boolean: Force redo the mean CNAs ans SNVs rates over individuals. Default FALSE.
#' @param forceRedoMatrixPlot boolean: Force redo the hit matrix plot. Default FALSE.
#' @param forceRedoMeanPlot boolean: Force redo the mean CNA plot. Default FALSE.
#' @param genome character: the genome being studied. Default "hg19".
#'
#' @details This function runs a cohort analysis, comparing two subgroups (within a project) to each other. See cohortAnalyseBatch for details.
#'
#'
compareGroups = function(metaDataFile, outputDirectories, project, subgroups1, subgroups2, name, includeNormal=F, clonalityCut=0.4, excludeSamples=c(), excludeIndividuals=c(), cosmicDirectory='', analysisName='cohortAnalysis', cpus=1, forceRedoVariants=F, forceRedoMean=F,
  forceRedoMeanPlot=F, forceRedoMatrixPlot=F, genome='hg19', ignoreCNAonly=F, cnvWeight=1) {

  metaData =
    makeMetaDataFromBatch(metaDataFile, outputDirectories, analysisName=analysisName,
                          excludeSamples=excludeSamples, excludeIndividuals=excludeIndividuals)
  createDirectories(metaData)
  linkBams(metaData)
  bringAnnotation(metaData, genome)
  
  samples = inProject(metaData, project, includeNormal=F, onlyDNA=T)
  individuals = metaData$samples[samples,]$INDIVIDUAL

  saveFile = paste0(metaData$project[project,]$Rdirectory, '/meanCNV.Rdata')
  if ( file.exists(saveFile) & !forceRedoMean ) {
    catLog('loading analysed cohort data...')
    load(saveFile)
    catLog('done.\n')
  }
  else {
    variants = getProjectVariants(metaData, project, cpus=cpus, forceRedo=forceRedoVariants)
    meanCNV = getMeanCNV(metaData, samples, variants, genome=genome, ignoreCNAonly=ignoreCNAonly)
    
    catLog('Saving meanCNV...')
    save(meanCNV, file=saveFile)
    catLog('done.\n')
  }

  meanCNV1 = meanCNV
  for ( subgroup in subgroups1 )
    meanCNV1 = splitMeanCNV(metaData, meanCNV1, project=project, subgroup=subgroup, includeNormal=includeNormal)$inGroup
  meanCNV2 = meanCNV
  for ( subgroup in subgroups2 )
    meanCNV2 = splitMeanCNV(metaData, meanCNV2, project=project, subgroup=subgroup, includeNormal=includeNormal)$inGroup
  #meanCNVboth = splitMeanCNV(metaData, meanCNV, project=project, subgroup=union(subgroups1, subgroups2), includeNormal=includeNormal)$inGroup
  meanCNVs = list('inGroup'=meanCNV1, 'outGroup'=meanCNV2)

  plotDirectory = paste0(metaData$project[project,]$plotDirectory, '/', name)
  ensureDirectoryExists(plotDirectory)

  plotFile = paste0(plotDirectory, '/meanCNVs.pdf')
  if ( !file.exists(plotFile) | forceRedoMeanPlot ) {
    pdf(plotFile, width=20, height=10)
    outputOut = plotMeanCNV(metaData, meanCNV=meanCNVs$outGroup, cosmicDirectory=cosmicDirectory,
      add=F, printGeneNames=F, meanCNV2=meanCNVs$inGroup, genome=genome)
    outputIn = plotMeanCNV(metaData, meanCNV=meanCNVs$inGroup, cosmicDirectory=cosmicDirectory,
      add=T, printGeneNames=T, genome=genome)
    legend('topright', c(subgroups1[1], subgroups1[1], subgroups1[1], subgroups2[1], subgroups2[1], subgroups2[1]), col=mcri(c('blue', 'red', 'green', 'cyan', 'orange', 'green')), lwd=c(2, 2, 2, 10, 10, 10), bg='white')
    dev.off()
  }
  
  plotMeanCNVtoFile(metaData, project, meanCNVs$inGroup, forceRedo=forceRedoMeanPlot,
                    plotFile=paste0(plotDirectory, '/meanCNV-', subgroups1[1], '.pdf'), genome=genome)
  plotMeanCNVtoFile(metaData, project, meanCNVs$outGroup, forceRedo=forceRedoMeanPlot,
                    plotFile=paste0(plotDirectory, '/meanCNV-', subgroups2[1], '.pdf'), genome=genome)
  
  plotSubgroupMutationMatrix(metaData, meanCNVs, project, name, cosmicDirectory=cosmicDirectory,
                             nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot)    
  plotMultipageMutationMatrix(metaData, meanCNVs$inGroup, project, cosmicDirectory=cosmicDirectory,
                              nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot,
                              plotFile=paste0(plotDirectory, '/mutationMatrix-', subgroups1[1], '.pdf'), genome=genome)
  plotMultipageMutationMatrix(metaData, meanCNVs$outGroup, project, cosmicDirectory=cosmicDirectory,
                              nGenes=30, pages=10, priorCnvWeight=cnvWeight, forceRedo=forceRedoMatrixPlot,
                              plotFile=paste0(plotDirectory, '/mutationMatrix-', subgroups2[1], '.pdf'), genome=genome)

  catLog('Done! Invisible returning mean CNV data.\n')
  invisible(meanCNV)
}


correctForGeneSize = function(hitMx, geneLengths) {

  
  
}
