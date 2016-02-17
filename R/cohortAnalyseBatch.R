

cohortAnalyseBatch =
  function(metaDataFile, outputDirectories, contrasts=list(), cpus=1, onlyDNA=T, clonalityCut=0.4,
           excludeSamples=c(), excludeIndividuals=c(), cosmicDirectory='',
           forceRedoVariants=F, forceRedoMean=F, forceRedoMatrixPlot=F, forceRedoMeanPlot=F, genome='hg19') {
  metaData =
    makeMetaDataFromBatch(metaDataFile, outputDirectories,
                          excludeSamples=excludeSamples, excludeIndividuals=excludeIndividuals)
  createDirectories(metaData)
  linkBams(metaData)
  bringAnnotation(metaData, genome)
  projects = getProjects(metaData, onlyDNA=onlyDNA)
  for ( project in projects ) {
    catLog('Cohort analysing', project, '\n')
    projectMeanCNV(metaData, project, cpus=cpus, onlyDNA=onlyDNA, clonalityCut=clonalityCut,
                   forceRedoMean=forceRedoMean, forceRedoVariants=forceRedoVariants, cosmicDirectory=cosmicDirectory,
                   forceRedoMatrixPlot=forceRedoMatrixPlot, forceRedoMeanPlot=forceRedoMeanPlot, genome=genome)
  }
}


makeMetaDataFromBatch =
  function(metaDataFile, outputDirectories, excludeSamples=c(), excludeIndividuals=c(),
           defaultProjectName='myProject', defaultCaptureRegions='myCapture',
           defaultBatchName='thisBatch', defaultDatatype='DNA') {
  samples = importSampleMetaData(metaDataFile)
  if ( any(samples$NAME %in% excludeSamples) ) {
    catLog('Excluding samples: ')
    catLog(samples[samples$NAME %in% excludeSamples,]$NAME, sep='\n')
    catLog('\n')
    samples = samples[!(samples$NAME %in% excludeSamples),]
  }
  if ( any(samples$INDIVIDUAL %in% excludeIndividuals) ) {
    catLog('Excluding individuals: ')
    catLog(unique(samples[samples$INDIVIDUAL %in% excludeIndividuals,]$INDIVIDUAL))
    catLog('\n')
    samples = samples[!(samples$INDIVIDUAL %in% excludeIndividuals),]
  }
  dataDirectory = paste0(dirname(metaDataFile), '/cohortAnalysis')
  if ( !file.exists(dataDirectory) ) dir.create(dataDirectory)
  
  samples$NORMAL = samples$NORMAL == 'YES'
  samples$RDIRECTORY = paste0(dataDirectory, '/R/samples/',samples$NAME)
  samples$PLOTDIRECTORY = paste0(dataDirectory, '/plots/samples/',samples$NAME)
  rownames(samples) = samples$NAME
  if ( !('PROJECT' %in% names(samples)) ) samples$PROJECT = rep(defaultProjectName, nrow(samples))
  if ( !('PROJECT.SUBGROUP' %in% names(samples)) ) samples$PROJECT.SUBGROUP = rep('', nrow(samples))
  if ( !('BATCH' %in% names(samples)) ) samples$BATCH = rep(defaultBatchName, nrow(samples))
  if ( !('DATATYPE' %in% names(samples)) ) samples$DATATYPE = rep(defaultDatatype, nrow(samples))
  if ( !('CAPTUREREGIONS' %in% names(samples)) ) samples$CAPTUREREGIONS = rep(defaultCaptureRegions, nrow(samples))

  
  individuals = data.frame('individual'=unique(samples$INDIVIDUAL))
  individuals$Rdirectory = paste0(dataDirectory, '/R/individuals/',individuals$NAME)
  individuals$plotDirectory = paste0(dataDirectory, '/plots/individuals/',individuals$NAME)
  rownames(individuals) = individuals$individual

  projects = if( 'PROJECT' %in% names(samples) ) {
      samples$PROJECT = make.names(samples$PROJECT)
      unique(unlist(strsplit(make.names(samples$PROJECT), split='\\.')))
    } else defaultProjectName
  projects = data.frame(
    'project'=projects,
    'Rdirectory'=paste0(dataDirectory, '/R/projects/',projects),
    'plotDirectory'=paste0(dataDirectory, '/plots/projects/',projects), stringsAsFactors=F)
  rownames(projects) = projects$project
  
  paths = list(
    'sources'=metaDataFile,
    'dataDirectory'=dataDirectory,
    'outputDirectories'=outputDirectories)

  ret = list('samples'=samples,
             'normalSamples'=samples[samples$NORMAL,],
             'individuals'=individuals,
             'project'=projects,
             'paths'=paths)
}


loadCohortMethods = function() {
  source('metaData.R')
  source('mergeMethods.R')
  source('meanCNV.R')
}


bringAnnotation = function(metaData, genome) {
  fromPath = paste0(metaData$paths$outputDirectories$Rdirectory, '/ensembl', genome, 'annotation.Rdata')
  toPath = paste0(metaData$paths$dataDirectory, '/resources/ensembl', genome, 'annotation.Rdata')
  if ( file.exists(fromPath) )
    system(paste0('cp ', fromPath, ' ', toPath))
}
