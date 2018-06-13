getCNV = function(metaData, samples) {
  Rdirs = unique(samplesToRdirs(samples, metaData))
  cnvs = list()
  for ( Rdir in Rdirs ) {
    catLog('Loading ', paste0(Rdir, '/clusters.Rdata'), '...', sep='')
    load(file=paste0(Rdir, '/clusters.Rdata'))
    catLog('done.\n')
    names = names(clusters)
    names = gsub('[_-]', '.', names)
    names(clusters) = gsub('[_-]', '.', names(clusters))
    names = names[names %in% samples]
    cnvs[names] = clusters[names]
  }
  cnvs = cnvs[samples]
  clusterColumns = Reduce(intersect, lapply(cnvs, function(cnv) colnames(cnv$clusters)))
  crColumns = Reduce(intersect, lapply(cnvs, function(cnv) colnames(cnv$CR)))
  cnvs = lapply(cnvs, function(cnv) {
    cnv$clusters =  cnv$clusters[,clusterColumns]
    cnv$CR =  cnv$CR[,crColumns]
    return(cnv)
  })
  return(cnvs)
}

getFit = function(metaData, samples) {
  Rdirs = unique(samplesToRdirs(samples, metaData))
  fits = list()
  for ( Rdir in Rdirs ) {
    loadFile = paste0(Rdir, '/fit.Rdata')
    if ( !file.exists(loadFile) ) loadFile = paste0(Rdir, '/fitP.Rdata')
    catLog('Loading ', loadFile, '...', sep='')
    fitName = load(file=loadFile)
    catLog('done.\n')
    if ( fitName == 'fit' ) fitP = fit$fit
    fitNames = gsub('-normal', '', colnames(fitP))
    standardNames = gsub('[_-]', '.', fitNames)
    fitNames = fitNames[standardNames %in% samples]
    standardNames = standardNames[standardNames %in% samples]
    fits[standardNames] = lapply(fitNames, function(name) subsetFit(fitP, cols=paste0(name, '-normal')))
  }
  fits = fits[samples]
  return(fits)
}

getVariant = function(metaData, samples, fillInMissing=F, cpus=1) {
  Rdirs = unique(samplesToRdirs(samples, metaData))
  variants = list('variants'=list(), 'SNPs'=c())
  for ( Rdir in Rdirs ) {
    a =
      try({
        catLog('Loading ', paste0(Rdir, '/allVariants.Rdata'), '...', sep='')
        load(file=paste0(Rdir, '/allVariants.Rdata'))
        catLog('done.\n')
        allVariants$variants$variants = cleanVariantRownames(allVariants$variants$variants)
        allVariants$variants$SNPs = allVariants$variants$SNPs[!duplicated(allVariants$variants$SNPs$x),]
        names = names(allVariants$variants$variants)
        names = gsub('[_-]', '.', names)
        names = names[names %in% samples]
        if ( length(names) == 0 ) stop('Couldnt find any relevant samples in ', paste0(Rdir, '/allVariants.Rdata'))
        names(allVariants$variants$variants) = gsub('[_-]', '.', names(allVariants$variants$variants))
        allVariants$variants$variants = allVariants$variants$variants[names]
        #don't bother with variants that arent supported by at least two reads
        presentVariants = unique(do.call(c, lapply(allVariants$variants$variants, function(q) rownames(q[q$var > 1,]))))
        allVariants$variants$variants = lapply(allVariants$variants$variants, function(q) q[presentVariants[presentVariants %in% rownames(q)],])
        if ( length(variants$variants) > 0 ) {
          if ( fillInMissing ) {
            variants = fillInMissingVariants(metaData, variants, allVariants$variants, cpus=cpus)
            allVariants$variants = fillInMissingVariants(metaData, allVariants$variants, variants, cpus=cpus)
          }
          allVariants$variants = matchVariants(allVariants$variants, variants)
          variants = matchVariants(variants, allVariants$variants)
          variants$variants[names] = allVariants$variants$variants
        }
        else variants = allVariants$variants
      })
    if ( class(a) == 'try-error' ) {
      warning('Failed to laod and merge variants from ', Rdir, '. Caught error message: ', a)

    }
  }
  samples = samples[samples %in% names(variants$variants)]
  variants$variants = variants$variants[samples]
  return(variants)
}

# This one doesnt work. May need to call missing variants in normals
# or retrieve from the saved variants in the normal folder.
getNormalVariant = function(metaData, variants, fillInMissing=F, cpus=1) {
  samples = names(variants$variants)
  Rdirs = unique(samplesToRdirs(samples))

  neededVariants = unique(do.call(c, lapply(variants$variants, rownames)))
  neededX = unique(variants$SNPs$x)
  ret = list(SNPs='', variants=list())
  for ( Rdir in Rdirs ) {
    #retrieve the normal variants
    load(paste0(Rdir, '/allVariants.Rdata'))
    normalVariants = allVariants$normalVariants

    #keep only needed variants
    normalVariants$variants = lapply(normalVariants$variants, function(q) q[rownames(q) %in% neededVariants,])
    normalVariants$SNPs = normalVariants$SNPs[normalVariants$SNPs$x %in% neededX,]

    ret$variants = c(ret$variants, normalVariants$variants)
    if ( class(ret$SNPs) == 'character' ) ret$SNPs = normalVariants$SNPs
    else {
      ret$SNPs = rbind(ret$SNPs, normalVariants$SNPs)
      ret$SNPs = ret$SNPs[!duplicated(ret$SNPs),]
    }
  }

  return(ret)
}

sampleToRdir = function(bamFile, metaData, RNA=F) {
  if ( 'outputDirectories' %in% names(metaData$paths) )
    return(metaData$paths$outputDirectories$Rdirectory)
  if ( gsub('^.+/', '', dirname(bamFile)) == 'bam' ) {
    if ( RNA ) {
      candidate = gsub('bam$', 'RNAcalling/R', dirname(bamFile))
      if ( file.exists(candidate) ) return(candidate)
    }
    candidate = gsub('bam$', 'R', dirname(bamFile))
    if ( file.exists(candidate) ) return(candidate)
  }
  candidate = paste0(dirname(bamFile), '/R')
  if ( RNA ) candidate = gsub('bam$', 'RNAcalling/R', dirname(bamFile))
  if ( file.exists(candidate) ) return(candidate)
  warning(paste0('Couldnt find R directory matching bam file', bamFile))
  return(NA)
}
samplesToRdirs = function(samples, metaData) {
  sapply(samples, function(sample) sampleToRdir(metaData$samples[sample,]$BAM, metaData,
                                                metaData$samples[sample,]$DATATYPE == 'RNA'))
}

updateScatterPlots = function(metaData, cosmicDirectory='', individuals=NA, forceRedo=F, forceRedoVariants=F, cpus=1, genome='hg19') {
  if ( is.na(individuals[1]) ) individuals = metaData$individuals$individual
  
  for ( individual in individuals ) {
    logFile = paste0(metaData$paths$dataDirectory, '/R/individuals/', individual, '/runtimeTracking.log')
    assign('catLog', function(...) {cat(..., file=logFile, append=T); cat(...)}, envir = .GlobalEnv)
    catLog('\n\n\n\n\n\n\n#############################################################\n')
    catLog('Running superFreq version', superVersion(), '\n')
    catLog('\n\nDoing', individual, '...\n')
    
    
    samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL == individual]
    normals = metaData$samples[samples,]$NORMAL
    names(normals) = samples

    Rdirectory = paste0(metaData$paths$dataDirectory, '/R/individuals/', individual)
    plotDirectory = paste0(metaData$paths$dataDirectory, '/plots/individuals/', individual)

        
    timePoints = metaData$samples[samples,]$TIMEPOINT
    names(timePoints) = samples
    dominantCaptureRegion = names(sort(table(metaData$samples$CAPTUREREGIONS[metaData$samples$INDIVIDUAL == individual]),
      decreasing=T))[1]
    genome = captureRegionsToGenome(dominantCaptureRegion)

    timeSeries = list(samples)
    names(timeSeries) = individual
    if ( length(timeSeries[[1]]) < 2 ) timeSeries = list()
    
    variants = getIndividualVariants(metaData, individual, cpus=cpus, forceRedo=forceRedoVariants)
    cnvs = getCNV(metaData, samples)
    stories = getStories(variants=variants, cnvs=cnvs, timeSeries=timeSeries,
      normals=normals, genome=genome, Rdirectory=Rdirectory, plotDirectory=plotDirectory, cpus=cpus, forceRedo=forceRedo)



    outputSomaticVariants(variants, genome=genome, plotDirectory=plotDirectory, cpus=cpus, forceRedo=forceRedo)
    variants = runVEP(variants, plotDirectory, cpus=cpus, genome=genome, forceRedoVEP=forceRedo)
    variants = getMoreVEPinfo(variants, plotDirectory, genome, cosmicDirectory=cosmicDirectory)

    samplePairs = metaToSamplePairs(samples, rep(individual, length(samples)), normals)
    inds = rep(individual, length(variants$variants))
    names(inds) = names(variants$variants)
    
    outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, forceRedo=forceRedo)
    makeSNPprogressionPlots(variants, timeSeries=timeSeries, normals = normals, plotDirectory=plotDirectory, forceRedo=forceRedo)
    makeRiverPlots(stories$stories, variants, genome=genome, cpus=cpus, plotDirectory=plotDirectory, forceRedo=forceRedo)
    makeScatterPlots(variants, samplePairs, timePoints, plotDirectory, genome=genome, cpus=cpus, forceRedo=forceRedo)
    makeCloneScatterPlots(variants, stories$stories, samplePairs, inds, timePoints,
                          plotDirectory, genome=genome, cpus=cpus, forceRedo=forceRedo)
    outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, forceRedo=forceRedo)


    
  }

  catLog('Done!\n')

}

updateCNVplots = function(metaData, forceRedo=F, cpus=1, genome='hg19') {
  samples = metaData$individuals$samples
  
  redoPairs = list()
  for ( ind in individuals ) {
    try({
      if ( sum(metaData$samples$INDIVIDUAL == ind) < 2 ) next
      samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL == ind]
      pairs = list()
      for ( i in 2:length(samples) )
        for ( j in 1:(i-1) )
          pairs = c(pairs, list(samples[c(j,i)]))
      
      for ( pair in pairs ) {
        freqDir = paste0(metaData$paths$dataDirectory, '/plots/samples/', pair[1], '/frequencyScatters')
        ensureDirectoryExists(freqDir)
        plotDir = paste0(freqDir, '/', pair[2])
        ensureDirectoryExists(plotDir)
        if ( forceRedo | length(list.files(plotDir)) == 0 ) redoPairs = c(redoPairs, list(pair))
      }
    })
  }

  if ( length(redoPairs) == 0 ) {
    catLog('Scatter plots were already in place.\n')
    return()
  }
  
  variants = getVariant(metaData, unique(unlist(redoPairs)), cpus=cpus)

  for ( pair in redoPairs ) {
    try({
      boring = (variants$variants[[pair[1]]]$var == 0 & variants$variants[[pair[2]]]$var == 0)
      q1 = variants$variants[[pair[1]]][!boring,]
      q2 = variants$variants[[pair[2]]][!boring,]
      catLog('Caculating flagged p-values...')
      ps=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, verbose=F)
      catLog(' and unflagged...')
      psuf=qualityScatter(q1, q2, variants$SNPs, cpus=cpus, plotFlagged=F, verbose=F)
      catLog('done.\n')
      
      freqDir = paste0(metaData$paths$dataDirectory, '/plots/samples/', pair[1], '/frequencyScatters')
      ensureDirectoryExists(freqDir)
      plotDir = paste0(freqDir, '/', pair[2])
      ensureDirectoryExists(plotDir)
      
      freqDir = paste0(metaData$paths$dataDirectory, '/plots/samples/', pair[2], '/frequencyScatters')
      ensureDirectoryExists(freqDir)
      linkDir = paste0(freqDir, '/', pair[1])
      ensureDirectoryExists(linkDir)
      
      outfile = paste0(plotDir, '/all.png')
      linkfile = paste0(linkDir, '/all.png')
      catLog('Plotting to', outfile, '\n')
      png(outfile, width = 10, height=10, res=144, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1]),
                     ylab=paste0('variant frequency for ', pair[2]), cpus=1)
      dev.off()
      if ( file.exists(linkfile) ) system(paste0('rm ', linkfile))
      system(paste0('ln -s ', outfile, ' ', linkfile))
      
      outfile = paste0(plotDir, '/allNamed.png')
      linkfile = paste0(linkDir, '/allNamed.png')
      catLog('Plotting to', outfile, '\n')
      png(outfile, width = 10, height=10, res=144, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=psuf, plotFlagged=F,
                     main=paste0('clean variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1]),
                     ylab=paste0('variant frequency for ', pair[2]), cpus=1,
                     print=T, printRedCut=0.25)
      dev.off()
      if ( file.exists(linkfile) ) system(paste0('rm ', linkfile))
      system(paste0('ln -s ', outfile, ' ', linkfile))
      
      outfile = paste0(plotDir, '/allFlagged.png')
      linkfile = paste0(linkDir, '/allFlagged.png')
      catLog('Plotting to', outfile, '\n')
      png(outfile, width = 10, height=10, res=144, units='in')
      qualityScatter(q1, q2, variants$SNPs, verbose=F, ps=ps,
                     main=paste0('all variants: ', pair[1], ' vs ', pair[2]),
                     xlab=paste0('variant frequency for ', pair[1]),
                     ylab=paste0('variant frequency for ', pair[2]), cpus=1)
      dev.off()
      if ( file.exists(linkfile) ) system(paste0('rm ', linkfile))
      system(paste0('ln -s ', outfile, ' ', linkfile))
      
      catLog('Now by chromsome. Preparing..')
      chrs = xToChr(q1$x, genome=genome)
      catLog('done!\n  Plotting chr:')
      for ( chr in names(chrLengths(genome)) ) {
        outfile = paste0(plotDir, '/chr', chr, '.png')
        linkfile = paste0(linkDir, '/chr', chr, '.png')
        png(outfile, width = 10, height=10, res=144, units='in')
        catLog(chr, '..', sep='')
        use = chrs == chr
        qualityScatter(q1[use,], q2[use,], variants$SNPs, ps=ps[use],
                       main=paste0('all variants: ', pair[1], ' vs ', pair[2], ', chr', chr),
                       xlab=paste0('variant frequency for ', pair[1]),
                       ylab=paste0('variant frequency for ', pair[2]), cpus=1, print=T,
                       printRedCut=0.25, plotPosition=T, verbose=F)
        dev.off()
        if ( file.exists(linkfile) ) system(paste0('rm ', linkfile))
        system(paste0('ln -s ', outfile, ' ', linkfile))
      }
      catLog('done!\n')
    })
  }

  for ( ind in individuals ) {
    indDir = paste0(metaData$paths$dataDirectory, '/plots/individuals/', ind)
    ensureDirectoryExists(indDir)
    outfile = paste0(indDir, '/newVariants.xls')
    
    samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL == ind]
    #set up links to sample plot directories
    for ( sample in samples ) {
      samplePlotDir = paste0(metaData$paths$dataDirectory, '/plots/samples/', sample)
      ensureDirectoryExists(samplePlotDir)
      linkfile = paste0(indDir, '/', sample)
      if ( file.exists(linkfile) ) system(paste0('rm ', linkfile))
      system(paste0('ln -s ', samplePlotDir, ' ', linkfile))
    }

    #ignore if not at least 2 samples
    if ( length(samples) < 2 ) next

    #set up pairs
    pairs = list()
    for ( i in 2:length(samples) )
      for ( j in 1:(i-1) )
        pairs = c(pairs, list(samples[c(j,i)]))

    #get the significantly different variants
    if ( (!file.exists(outfile) | forceRedo) & length(pairs) > 0 ) {
      news = list()
      for ( pair in pairs ) {
        name1 = substring(gsub('\\.', '', paste0(pair[1], ' to ', pair[2])), 1, 31)
        name2 = substring(gsub('\\.', '', paste0(pair[2], ' to ', pair[1])), 1, 31)
        catLog('Looking for new cancer variants in ', name1, '\n')
        news[[name1]] = newVariants(variants$variants[[pair[1]]], variants$variants[[pair[2]]], variants$SNPs, genome, cpus=cpus)
        catLog('Looking for new cancer variants in ', name2, '\n')
        news[[name2]] = newVariants(variants$variants[[pair[2]]], variants$variants[[pair[1]]], variants$SNPs, genome, cpus=cpus)
      }

      #print to excel
      catLog('Writing to', outfile, '\n')
      WriteXLS('news', outfile)
    }
  }
  catLog('Done!\n')

}


updateMultisamplePlots = function(metaData, forceRedo=F) {
  individuals = metaData$individuals$individual
  
  redoInd = c()
  for ( ind in individuals ) {
    if ( sum(metaData$samples$INDIVIDUAL == ind) < 2 ) next
    indDir = paste0(metaData$paths$dataDirectory, '/plots/individuals/', ind)
    ensureDirectoryExists(indDir)
    plotFile = paste0(indDir, '/multisample.pdf')
    if ( forceRedo | !file.exists(plotFile) ) redoInd = c(redoInd, ind)
  }
  
  if ( length(redoInd) == 0 ) {
    catLog('Multisample plots were already in place.\n')
    return()
  }
  
  variants = getVariant(metaData, metaData$samples$NAME[metaData$samples$INDIVIDUAL %in% redoInd])
  
  for ( ind in redoInd ) {
    try({
      indDir = paste0(metaData$paths$dataDirectory, '/plots/individuals/', ind)
      timeSeries = metaData$samples$NAME[metaData$samples$INDIVIDUAL == ind]
      normals = metaData$samples$NORMAL[metaData$samples$INDIVIDUAL == ind]
      if ( length(timeSeries) < 2 ) next
      outfile = paste0(indDir, '/frequencyHeatmap.pdf')
      excelFileDB = paste0(indDir, '/frequencyHeatmap_DB.xls')
      excelFileNotDB = paste0(indDir, '/frequencyHeatmap_NotDB.xls')
      catLog('Plotting SNP progression to ', outfile, '..', sep='')
      pdf(outfile, width = 15, height=10)
      qualityProgression(variants$variants[timeSeries], variants$SNPs, normals, nondb=F,
                         excelFile=excelFileDB, main='dbSNPs only')
      qualityProgression(variants$variants[timeSeries], variants$SNPs, normals, db=F,
                         excelFile=excelFileNotDB, main='non-dbSNPs only')
      dev.off()
      catLog('done.\n')
    })
  }
  catLog('Done!\n')

}

cleanVariants = function(variants) {
  variants$variants = lapply(variants$variants, function(q) {
    q = q[!is.na(q$x),]
    q$severity[is.na(q$severity)] = 100
    q$type[is.na(q$type)] = 'notChecked'
    return(q)
    })
  return(variants)
}


getProjectVariants = function(metaData, project, cpus=1, onlyDNA=T, includeNormal=F, forceRedo=F) {
  if ( length(inProject(metaData, project, includeNormal=T, onlyDNA=onlyDNA)) == 0 ) {
    warning(paste0('couldnt find any samples in project ', project, '.'))
    return()
  }
  saveFile = paste0(metaData$project[project,]$Rdirectory, '/variants.Rdata')
  ensureDirectoryExists(metaData$project[project,]$Rdirectory)
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading saved project variants...')
    load(saveFile)
    return(variants)
  }

  samples = inProject(metaData, project, includeNormal=includeNormal, onlyDNA=onlyDNA)
  variants = getVariant(metaData, samples, cpus=cpus)

  catLog('Saving project variants to', saveFile, '.\n')
  save(variants, file=saveFile)
  return(variants)
}

getAllIndividualVariants = function(metaData, individual, cpus=1, forceRedo=F) {
  if ( !(individual %in% metaData$individuals$individual) ) {
    warning(paste0('couldnt find individual ', individual, ' in metadata.'))
    return()
  }
  saveFile = paste0(metaData$individuals[individual,]$Rdirectory, '/allVariants.Rdata')
  ensureDirectoryExists(metaData$individuals[individual,]$Rdirectory)
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading saved variants...')
    load(saveFile)
    return(allVariants)
  }

  #import the variants from their respective batch directories. Match and cross-check variants if needed.
  variants = getIndividualVariants(metaData, individual, cpus=cpus, forceRedo=forceRedo)

  #get the variants in the individuals from the normal samples
  dominantCaptureRegion = names(sort(table(metaData$samples$CAPTUREREGIONS[metaData$samples$INDIVIDUAL == individual]),decreasing=T))[1]
  genome = captureRegionsToGenome(dominantCaptureRegion)
  BQoffset = captureRegionsToBQoffset(dominantCaptureRegion)
  Rdirectory = paste0(metaData$paths$dataDirectory, '/R/individuals/', individual)
  plotDirectory = paste0(metaData$paths$dataDirectory, '/plots/individuals/', individual)
  
  normalPath = captureRegionsToNormalPath(dominantCaptureRegion)
  normalBamFiles = list.files(paste0(normalPath, '/bam'), pattern = '*.bam$', full.names=T)
  normalNames = gsub('^.*/', '', gsub('.bam$', '', normalBamFiles))
  normalRdirectory = paste0(normalPath, '/R')
  normalVariants = getNormalVariants(variants, normalBamFiles, normalNames, dominantCaptureRegion,
                                     genome, BQoffset, normalRdirectory, Rdirectory, plotDirectory, cpus=cpus,
                                     forceRedoSNPs=F, forceRedoVariants=F)

  #redo flagging from the normals, and somatic analysis, by first resetting whatever was done in the batch
  #In particular this will allow the somatic score to be calculated using a matched normal from another batch.
  flagsFromNormal = c('Nnn', 'Nnc', 'Vn', 'Mc', 'Pn', 'Nr')
  variants$variants = lapply(variants$variants, function(q) {
    for ( flag in flagsFromNormal )
      q$flag = gsub(flag, '', q$flag)
    q$somaticP = 0
    return(q)
  })

  samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL == individual]
  normals = metaData$samples[samples,]$NORMAL
  names(normals) = samples
  allVariants = matchFlagVariants(variants, normalVariants, individual, normals, Rdirectory,
                                  cpus=cpus, forceRedoMatchFlag=forceRedo)

  return(allVariants)
}


getIndividualVariants = function(metaData, individual, cpus=1, forceRedo=F) {
  if ( !(individual %in% metaData$individuals$individual) ) {
    warning(paste0('couldnt find individual ', individual, ' in metadata.'))
    return()
  }
  saveFile = paste0(metaData$individuals[individual,]$Rdirectory, '/variants.Rdata')
  ensureDirectoryExists(metaData$individuals[individual,]$Rdirectory)
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading saved variants...')
    load(saveFile)
    return(variants)
  }

  samples = rownames(metaData$samples)[metaData$samples$INDIVIDUAL == individual]
  variants = getVariant(metaData, samples, fillInMissing=T, cpus=cpus)
  keepX = unique(do.call(c, lapply(variants$variants, function(q) q$x[q$var > 1])))
  variants$SNPs = variants$SNPs[variants$SNPs$x %in% keepX,]
  variants$variants = lapply(variants$variants, function(q) q[q$x %in% keepX,])
  keepVar = unique(do.call(c, lapply(variants$variants, function(q) rownames(q)[q$var > 1])))
  variants$variants = lapply(variants$variants, function(q) q[keepVar,])
    
  save(variants, file=saveFile)

  return(variants)
}


getIndividualSmallVariants = function(metaData, individual, cpus=1, forceRedo=F) {
  if ( !(individual %in% metaData$individuals$individual) ) {
    warning(paste0('couldnt find individual ', individual, ' in metadata.'))
    return()
  }

  saveFile = paste0(metaData$individuals[individual,]$Rdirectory, '/smallVariants.Rdata')
  ensureDirectoryExists(metaData$individuals[individual,]$Rdirectory)
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading saved small variants...')
    load(saveFile)
    return(variants)
  }

  variants =
    getIndividualVariants(metaData=metaData, individual=individual,
                          cpus=cpus, forceRedo=forceRedo)

  qs = variants$variants
  keep = sapply(qs, function(q) q$flag=='' & q$cov > 5 & q$var/q$cov > 0.05 & !is.na(q$cov))
  keep = apply(keep, 1, any)

  keepColumns = c('x', 'cov', 'ref', 'var', 'flag', 'db', 'dbMAF', 'dbValidated', 'somaticP', 'germline', 'severity', 'isCosmicCensus')
  smallVariants = lapply(qs, function(q) q[keep,names(q) %in% keepColumns])

  save(smallVariants, file=saveFile)
  return(smallVariants)
}

mergeStrs = function(strs, sep='') {
  return(do.call(paste, c(as.list(strs), sep='')))
}

#takes all the variants present in qs2 that are not present in qs1,
#checks them up in the BAMs of qs1 and adds them to qs1, which is returned.
fillInMissingVariants = function(metaData, qs1, qs2, cpus=1) {
  #check that qs1 and qs2 share columns.
  if ( !setequal(colnames(qs1$variants[[1]]), colnames(qs2$variants[[1]])) ) {
    warning(paste('FillInMissingVariants called with variants with different columns. qs1 has', mergeStrs(colnames(qs1$variants[[1]])), ', qs2 has', mergeStrs(colnames(qs2$variants[[1]])), '.'))
    commonCols = intersect(colnames(qs1$variants[[1]]), colnames(qs2$variants[[1]]))
    qs1$variants = lapply(qs1$variants, function(q) q[,commonCols])
    qs2$variants = lapply(qs2$variants, function(q) q[,commonCols])
  }
  
  #if columns not in the same order, have qs1 match the qs2 order.
  if ( !identical(colnames(qs1$variants[[1]]), colnames(qs2$variants[[1]])) ) {
    catLog('Matching order of columns.\n')
    qs1$variants = lapply(qs1$variants, function(q) q[,colnames(qs2$variants[[1]],)])
  }
  
  allSNPs = qs2$SNPs[!duplicated(qs2$SNPs$x),]
  q2 = qs2$variants[[1]]
  q2 = q2[!duplicated(paste0(q2$x, q2$variant)),]
  q2 = q2[order(q2$x, q2$variant),]
  allX = q2$x
  allVar = q2$variant
  allRownames = rownames(q2)

  catLog('Filling in missing variants: Got', length(allRownames), 'variants that need to be present.\n')

  for ( i in 1:length(qs1$variants) ) {
    #first find the missing variants
    name = names(qs1$variants)[i]
    q = qs1$variants[[i]]
    bam = metaData$samples[names(qs1$variants)[i],]$BAM
    missing = !(allRownames %in% rownames(q))
    catLog('Sample', names(qs1)[i], 'has', sum(missing), 'missing variants.\n')
    if ( sum(missing) == 0 ) next

    #see if the positions of the missing variants are already called in q1 (for other variants)
    #then we can re-use the coverage and set the variant count to 0, dont need to re-check BAM
    positionChecked = allX[missing] %in% q$x
    catLog(sum(positionChecked), 'can be reused from previous analysis.\n')
    if ( sum(positionChecked) > 0 )
      q = rbind(q, shareVariants(list(q[q$x %in% allX[missing][positionChecked],], q2[missing,][positionChecked,]))[[1]])
    


    #for the remaining variants, we need to check the BAM for q1.
    missing = !(allRownames %in% rownames(q))
    catLog(sum(missing), 'need to be called from bam.\n')
    if ( sum(missing) > 0 ) {
      newSNPs = allSNPs[allSNPs$x %in% allX[missing],]
      BQoffset = captureRegionsToBQoffset(metaData$samples[name,]$CAPTUREREGIONS)
      genome = captureRegionsToGenome(metaData$samples[name,]$CAPTUREREGIONS)
      newQ = bamToVariants(bam=bam, positions=newSNPs, BQoffset=BQoffset, genome=genome, cpus=cpus)[[1]]
      #newQ = QCsnps(pileups=importQualityScores(newSNPs, bam, BQoffset, genome=genome, cpus=cpus)[[1]],
      #  SNPs=newSNPs, cpus=cpus)
      catLog('Matching new calls to the missing variants.\n')
      newQ = newQ[order(newQ$x, newQ$variant),]
      #fill in empty entry if the right variant wasn't called over the position
      newQ = rbind(newQ, shareVariants(list(newQ, q2[missing,]))[[1]])
      #only keep the missing variants that are present in q2 but not q1
      newQ = newQ[rownames(newQ) %in% allRownames[missing],]

      #add any missing columns, such as db, VEP info, etc.
      missingColumns = colnames(q2)[!(colnames(q2) %in% colnames(newQ))]
      if ( length(missingColumns) > 0 ) {
        newQ = cbind(newQ, q2[rownames(newQ),missingColumns])
      }
      #excpet somaticP that should be set to 0, as it wasnt called in q1.
      if ( 'somaticP' %in% colnames(newQ) ) newQ$somaticP = 0

      #any missing information, such as gene, dbSNP information or VEP output can be copied from q2.
      missingColumns = colnames(q2)[!(colnames(q2) %in% newQ)]
      if ( length(missingColumns) > 0 )
        newQ = cbind(newQ, q2[rownames(newQ),missingColumns])
      newQ = newQ[,colnames(q)]
      
      catLog('Adding', nrow(newQ), 'freshly called variants.\n')
      q = rbind(q, newQ)
    }
    
    q = q[order(q$x, q$variant),]
    qs1$variants[[i]] = q
  }

  if ( length(qs1$SNPs) != length(qs2$SNPs) || any(colnames(qs2$SNPs) != colnames(qs1$SNPs)) ) {
    catLog('Matching SNP columns.\n')
    commonColumns = intersect(colnames(qs1$SNPs), colnames(qs2$SNPs))
    qs1$SNPs = qs1$SNPs[,commonColumns]
    qs2$SNPs = qs2$SNPs[,commonColumns]
  }
  
  
  if ( any(!(qs2$SNPs$x %in% qs1$SNPs$x)) ) {
    catLog('Matching meta-information about variants.\n')
    qs1$SNPs = rbind(qs1$SNPs, qs2$SNPs[!(qs2$SNPs$x %in% qs1$SNPs$x),])
  }
  
  return(qs1)
}


updateStory = function(metaData, individual, cpus=1, forceRedo=F) {
  samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL==individual]
  timeSeries = list(samples)
  names(timeSeries) = individual
  normals = metaData$samples[samples,]$NORMAL
  names(normals) = samples

  variants = getIndividualVariants(metaData, individual)
  dominantCaptureRegion = names(sort(table(metaData$samples$CAPTUREREGIONS[metaData$samples$INDIVIDUAL == individual]),decreasing=T))[1]
  genome = captureRegionsToGenome(dominantCaptureRegion)
  BQoffset = captureRegionsToBQoffset(dominantCaptureRegion)
  Rdirectory = paste0(metaData$paths$dataDirectory, '/R/individuals/', individual)
  plotDirectory = paste0(metaData$paths$dataDirectory, '/plots/individuals/', individual)
  
  normalPath = captureRegionsToNormalPath(dominantCaptureRegion)
  normalBamFiles = list.files(paste0(normalPath, '/bam'), pattern = '*.bam$', full.names=T)
  normalNames = gsub('^.*/', '', gsub('.bam$', '', normalBamFiles))
  normalRdirectory = paste0(normalPath, '/R')
  
  normalVariants = getNormalVariants(variants, normalBamFiles, normalNames, dominantCaptureRegion,
                                     genome, BQoffset, normalRdirectory, Rdirectory, plotDirectory, cpus=cpus,
                                     forceRedoSNPs=F, forceRedoVariants=F)
  
  cnvs = getCNV(metaData, samples)

  stories = getStories(variants, cnvs, timeSeries, normals, Rdirectory, plotDirectory, cpus=cpus, forceRedo=F)

  variants = getAllVEPdata(metaData, variants)

  
}

captureRegionsToNormalPath = function(captureRegions, coverage=F) {
  warning('Called a hardcoding function!')
  if ( captureRegions == 'SureSelect' & !coverage ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/Agilent')
  if ( captureRegions == 'SureSelect' & coverage ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/AgilentBlood')
  if ( captureRegions == 'TrueSeq' ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/Illumina')
  if ( captureRegions == 'RNA' ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/RNA')
  if ( captureRegions == 'kinome' ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/kinome')
  if ( captureRegions == 'genome' ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/genome')
  if ( captureRegions == 'siblingMUD' ) return('/wehisan/general/academic/grp_leukemia_genomics/normals/siblingMUD')
  warning(paste0('Dont know a normal path for capture regions ', captureRegions))
  return('')
}
captureRegionsToBQoffset = function(captureRegions) {
  warning('Called a hardcoding function!')
  if ( captureRegions == 'SureSelect' ) return(33)
  if ( captureRegions == 'TrueSeq' ) return(64)
  if ( captureRegions == 'RNA' ) return(33)
  if ( captureRegions == 'kinome' ) return(33)
  if ( captureRegions == 'genome' ) return(33)
  if ( captureRegions == 'siblingMUD' ) return(33)
  warning(paste0('Dont know a BQoffset for capture regions ', captureRegions))
  return('')
}
captureRegionsToGenome = function(captureRegions) {
  warning('Called a hardcoding function!')
  if ( captureRegions == 'SureSelect' ) return('hg19')
  if ( captureRegions == 'TrueSeq' ) return('hg19')
  if ( captureRegions == 'RNA' ) return('hg19')
  if ( captureRegions == 'kinome' ) return('hg19')
  if ( captureRegions == 'genome' ) return('hg19')
  if ( captureRegions == 'siblingMUD' ) return('hg19')
  warning(paste0('Dont know a genome for capture regions ', captureRegions))
  return('')
}
captureRegionsToBed = function(captureRegions) {
  warning('Called a hardcoding function!')
  if ( captureRegions == 'SureSelect' ) return('/wehisan/general/academic/grp_leukemia_genomics/AGILENT_EXOME/HUMAN_v5/hg19.bed')
  if ( captureRegions == 'TrueSeq' ) return('/wehisan/general/academic/grp_leukemia_genomics/F13TSFAPHT0461_HUMwkbX/analysis/captureRegions.bed')
  if ( captureRegions == 'RNA' ) return('/wehisan/general/academic/grp_leukemia_genomics/data/captureRegions/RNA/captureRegions.bed')
  if ( captureRegions == 'kinome' ) return('/wehisan/general/academic/grp_leukemia_genomics/DLBCL\\ PROJECT/analysis/captureRegions.gc.bed')
  if ( captureRegions == 'genome' ) return('/wehisan/general/academic/grp_leukemia_genomics/data/captureRegions/genomic/captureRegions.hg19.10k.bed')
  if ( captureRegions == 'siblingMUD' ) return('/wehisan/general/academic/grp_leukemia_genomics/data/captureRegions/siblingMUD/amplicons.bed')
  warning(paste0('Dont know a bed file path for capture regions ', captureRegions))
  return('')
}

getAllVEPdata = function(metaData, variants) {
  samples = names(variants$variants)
  plotDirectories = gsub('R$', 'plots', samplesToRdirs(samples, metaData))
  for ( sample in samples ) {
    variantSample = variants
    variantSample$variants = variantSample$variants[sample]
    variantSample = getMoreVEPinfo(variantSample, plotDirectories[sample])
    variants$variants[[sample]] = variantSample$variants[[sample]]
  }
  return(variants)
}

#Returns the samples that are part of a project
getProjects = function(metaData, onlyDNA=T) {
  projects = strsplit(metaData$samples$PROJECT, ',')
  ret = unique(unlist(projects))
  return(ret)
}

#Returns the samples that are part of a project
inProject = function(metaData, project, includeNormal=T, onlyDNA=T) {
  projects = strsplit(metaData$samples$PROJECT, ',')
  isInProject = sapply(projects, function(sampleProjects) any(sampleProjects == project))
  if ( !includeNormal ) isInProject = isInProject & !metaData$samples$NORMAL
  if ( onlyDNA ) isInProject = isInProject & metaData$samples$DATATYPE == 'DNA'
  return(metaData$samples$NAME[isInProject])
}

#This function returns the unique subgroups of the provided project
getSubgroups = function(metaData, project, includeNormal=T, onlyDNA=T) {
  samples = inProject(metaData, project, includeNormal=includeNormal, onlyDNA=onlyDNA)
  subgroups = metaData$samples[samples,]$PROJECT.SUBGROUP
  subgroups = subgroups[subgroups != '']
  subgroupList = strsplit(subgroups, ',')
  uniqueSubgroups = unique(unlist(subgroupList))
  isInformative = sapply(uniqueSubgroups, function(subgroup) !all(sapply(subgroupList, function(groups) subgroup %in% groups)) & any(sapply(subgroupList, function(groups) subgroup %in% groups)))
  return(uniqueSubgroups[isInformative])
}

#Returns the samples that are part of a project and subgroups
inSubgroup = function(metaData, project, subgroup, includeNormal=T, onlyDNA=T) {
  projectSamples = inProject(metaData, project, includeNormal=includeNormal, onlyDNA=onlyDNA)
  subgroups = metaData$samples[projectSamples,]$PROJECT.SUBGROUP
  subgroups = strsplit(subgroups, ',')
  isInSubgroup = sapply(subgroups, function(groups) any(groups %in% subgroup))
  subgroupSamples = projectSamples[isInSubgroup]
  if ( !includeNormal ) subgroupSamples = subgroupSamples[!metaData$samples[subgroupSamples,]$NORMAL]
  if ( onlyDNA ) subgroupSamples = subgroupSamples[metaData$samples[subgroupSamples,]$DATATYPE == 'DNA']
  return(subgroupSamples)
}

#returns the individual that the samples come from.
sampleToIndividual = function(metaData, samples) {
  return(metaData$samples[samples,]$INDIVIDUAL)
}

#takes a list of variants, and makes sure the rowname of each list entry is the appropriate combination
#of the position and variant. Returns the list of variants with fixed rownames.
cleanVariantRownames = function(variants) {
  if ( length(variants) == 0 ) return(variants)

  options(scipen = 10)
  variants = lapply(variants, function(q) {
    q = q[!is.na(q$x),]
    rownames(q) = paste0(q$x, q$variant)
    return(q)
    })
  return(variants)
}
