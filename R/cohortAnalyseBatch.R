
#' Analyse individuals for reccuring mutations
#'
#' @param metaDataFile character: path to the metaData file.
#' @param outputDirectories A named list of output directories, containing the entries Rdirectory and plotDirectory where the saved data and plots will be stored respectively.
#' @param cpus integer: the maximum number of cpus to run on.
#' @param clonalityCut numeric: the minimum required clonality to be included in the analysis. Deafult 0.4.
#' @param excludeSamples character: The samples to be excluded from the analysis. Default c().
#' @param excludeIndividuals character: The individuals to be excluded from the analysis. Default c().
#' @param cosmicDirectory character: The directory with the COSMIC data.
#' @param analysisName character: The name of the directory where the analysis results are saved.
#'                                Default "cohortAnalysis".
#' @param forceRedoVariants boolean: Force redo the merging of the variants between individuals. Default FALSE.
#' @param forceRedoMean boolean: Force redo the mean CNAs ans SNVs rates over individuals. Default FALSE.
#' @param forceRedoMatrixPlot boolean: Force redo the hit matrix plot. Default FALSE.
#' @param forceRedoMeanPlot boolean: Force redo the mean CNA plot. Default FALSE.
#' @param genome character: the genome being studied. Default "hg19".
#'
#' @details This function calculates mutation rates over genes, both protein changing SNVs, as well as CNA rates for complete loss, loss, gain (3 copies) and amplification (4 or more copies). It also track biallelic loss of genes in samples, by complete loss, a protein changing SNV plus a loss, or two SNVs (that are assumed to be on different alleles). Output is a plot over the genome of the CNA rates, as well as a "top table" of frequently mutated genes. It is run as an afterburner, and needs a finished superFreq analysis to be run on the samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metaDataFile = '/absolute/path/to/metaData.txt'
#'
#' Rdirectory = '/absolute/path/to/R'
#' plotDirectory = '/absolute/path/to/plots'
#'
#' cpus=6
#' genome = 'hg19'
#'
#' outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
#' 
#' cohortAnalyseBatch(metaDataFile, outputDirectories, cpus=cpus, genome=genome)
#'
#' }
cohortAnalyseBatch = function(metaDataFile, outputDirectories, cpus=1, onlyDNA=T, clonalityCut=0.4, includeNormal=F,
  excludeSamples=c(), excludeIndividuals=c(), cosmicDirectory='', analysisName='cohortAnalysis', cnvWeight=1,
  forceRedoVariants=F, forceRedoMean=F, forceRedoMatrixPlot=F, forceRedoMeanPlot=F, genome='hg19', ignoreCNAonly=F) {

  ensureDirectoryExists(outputDirectories$Rdirectory, verbose=F)
  logFile = paste0(normalizePath(outputDirectories$Rdirectory), '/runtimeTracking.log')
  assign('catLog', function(...) {cat(..., file=logFile, append=T); cat(...)}, envir = .GlobalEnv)
  catLog('Running superFreq version', superVersion(), '\n')

  metaData =
    makeMetaDataFromBatch(metaDataFile, outputDirectories, analysisName=analysisName,
                          excludeSamples=excludeSamples, excludeIndividuals=excludeIndividuals)
  createDirectories(metaData)
  linkBams(metaData)
  bringAnnotation(metaData, genome)
  projects = getProjects(metaData, onlyDNA=onlyDNA)
  for ( project in projects ) {
    catLog('Cohort analysing', project, '\n')
    a = projectMeanCNV(metaData, project, cpus=cpus, onlyDNA=onlyDNA, clonalityCut=clonalityCut, includeNormal=includeNormal,
      forceRedoMean=forceRedoMean, forceRedoVariants=forceRedoVariants, cosmicDirectory=cosmicDirectory, cnvWeight=cnvWeight,
      forceRedoMatrixPlot=forceRedoMatrixPlot, forceRedoMeanPlot=forceRedoMeanPlot, genome=genome, ignoreCNAonly=ignoreCNAonly)
    if ( class(a) == 'try-error' ) {
      catLog('Failed project mean CNV for project ', project, ' with error message: ', a, '\n', sep='')
      warning('Failed project mean CNV for project ', project, ' with error message: ', a)
    }
  }

  invisible(metaData)
}


#' Analyse individuals for reccuring mutations
#'
#' @param metaDataFile character: path to the metaData file.
#' @param outputDirectories A named list of output directories, containing the entries Rdirectory and plotDirectory where the saved data and plots will be stored respectively.
#' @param project character: The project containing the subgroups.
#' @param subgroup1 character: The first subgroup(s).
#' @param subgroup2 character: The second subgroup(s).
#' @param name character: The name of the comparison. This names the output directory.
#' @param cpus integer: the maximum number of cpus to run on.
#' @param clonalityCut numeric: the minimum required clonality to be included in the analysis. Deafult 0.4.
#' @param excludeSamples character: The samples to be excluded from the analysis. Default c().
#' @param excludeIndividuals character: The individuals to be excluded from the analysis. Default c().
#' @param cosmicDirectory character: The directory with the COSMIC data.
#' @param analysisName character: The name of the directory where the analysis results are saved.
#'                                Default "cohortAnalysis".
#' @param forceRedoVariants boolean: Force redo the merging of the variants between individuals. Default FALSE.
#' @param forceRedoMean boolean: Force redo the mean CNAs ans SNVs rates over individuals. Default FALSE.
#' @param forceRedoMatrixPlot boolean: Force redo the hit matrix plot. Default FALSE.
#' @param forceRedoMeanPlot boolean: Force redo the mean CNA plot. Default FALSE.
#' @param genome character: the genome being studied. Default "hg19".
#'
#' @details This function calculates mutation rates over genes, both protein changing SNVs, as well as CNA rates for complete loss, loss, gain (3 copies) and amplification (4 or more copies). It also track biallelic loss of genes in samples, by complete loss, a protein changing SNV plus a loss, or two SNVs (that are assumed to be on different alleles). Output is a plot over the genome of the CNA rates, as well as a "top table" of frequently mutated genes. It is run as an afterburner, and needs a finished superFreq analysis to be run on the samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metaDataFile = '/absolute/path/to/metaData.txt'
#'
#' Rdirectory = '/absolute/path/to/R'
#' plotDirectory = '/absolute/path/to/plots'
#'
#' cpus=6
#' genome = 'hg19'
#'
#' outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
#' 
#' cohortAnalyseBatch(metaDataFile, outputDirectories, cpus=cpus, genome=genome)
#'
#' }
cohortAnalyseBatchContrast = function(metaDataFile, outputDirectories, project, subgroups1, subgroups2, name, cpus=1,
  onlyDNA=T, clonalityCut=0.4, 
  excludeSamples=c(), excludeIndividuals=c(), cosmicDirectory='', analysisName='cohortAnalysis', cnvWeight=1,
  forceRedoVariants=F, forceRedoMean=F, forceRedoMatrixPlot=F, forceRedoMeanPlot=F, genome='hg19', ignoreCNAonly=F) {

  logFile = normalizePath(paste0(Rdirectory, '/runtimeTracking.log'))
  assign('catLog', function(...) {cat(..., file=logFile, append=T); cat(...)}, envir = .GlobalEnv)
  catLog('Running superFreq version', superVersion(), '\n')

  metaData =
    compareGroups(metaDataFile=metaDataFile, outputDirectories=outputDirectories, project=project,
                  subgroups1=subgroups1, subgroups2=subgroups2,
                  name=name, clonalityCut=clonalityCut, excludeSamples=excludeSamples, excludeIndividuals=excludeIndividuals,
                  cosmicDirectory=cosmicDirectory, analysisName=analysisName, cpus=cpus, forceRedoVariants=forceRedoVariants,
                  forceRedoMean=forceRedoMean, ignoreCNAonly=ignoreCNAonly, cnvWeight=cnvWeight,
                  forceRedoMeanPlot=forceRedoMeanPlot, forceRedoMatrixPlot=forceRedoMatrixPlot, genome=genome)

  invisible(metaData)
}


makeMetaDataFromBatch =
  function(metaDataFile, outputDirectories, analysisName, excludeSamples=c(), excludeIndividuals=c(),
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
  dataDirectory = paste0(dirname(metaDataFile), '/', analysisName)
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
    unique(unlist(strsplit(samples$PROJECT, split=',')))
  } else defaultProjectName
  projects = data.frame(
    'project'=projects,
    'Rdirectory'=paste0(dataDirectory, '/R/projects/',projects),
    'plotDirectory'=paste0(dataDirectory, '/plots/projects/',projects), stringsAsFactors=F)
  rownames(projects) = projects$project
  
  paths = list(
    'sources'=normalizePath(metaDataFile),
    'dataDirectory'=normalizePath(dataDirectory),
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


#' Merges data from several batches for cohort analysis
#'
#' @param paths data.frame: One row for each batch to be merged.
#"                          Columns are metaDataFile and Rdirectory.
#' @param targetMetaDataFile character: Path to where the merged metaData will be stored.
#'                           This also affect where a downstream cohort analysis is stored.
#' @param targetRdirectory character: Path to where the Rdata will be stored.
#'
#' @details This function merges the data from multiple batches, setting them up for a joint cohort analysis. Note that this is only for the purpose of a downstream cohort analysis. The output R directory cannot be run as a batch Rdirectory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths = data.frame('metaDataFile' = c('~/superFreq/myFirstBatch/metaData.txt', '~/superFreq/mySecondBatch/metaData.txt'),
#'            'Rdirectory'= c('~/superFreq/myFirstBatch/R', '~/superFreq/myFirstBatch/R'))
#'
#' targetMetaDataFile = '~/superFreq/mergedAnalysis/metaData.txt'
#' targetRdirectory = '~/superFreq/mergedAnalysis/R'
#'
#' mergeBatches(paths, targetMetaDataFile, targetRdirectory)
#'
#' outputDirectories = list('Rdirectory'=targetRdirectory)
#' cohortAnalyseBatch(targetMetaDataFile, outputDirectories, cpus=6, genome='hg19')
#'
#' }
mergeBatches = function(paths, targetMetaDataFile, targetRdirectory, cpus=1) {  
  if ( !('metaDataFile' %in% names(paths)) ) stop('paths need to have an metaDataFile column.')
  if ( !('Rdirectory' %in% names(paths)) ) stop('paths need to have an Rdirectory column.')

  assign('catLog', function(...) {cat(...)}, envir = .GlobalEnv)
  catLog('Running superFreq version', superVersion(), '\n')

  variantsSaveFile = paste0(targetRdirectory, '/allVariants.Rdata')
  clustersSaveFile = paste0(targetRdirectory, '/clusters.Rdata')
  if ( file.exists(variantsSaveFile) & file.exists(clustersSaveFile) ) {
    catLog('Batches already merged. Moving along.\n')
    return()
  }

  
  #merge metaData files.
  metaDatas = lapply(as.character(paths$metaDataFile), function(metaDataFile) importSampleMetaData(metaDataFile))
  #fill in missing columns with blanks.
  columns = Reduce(union, lapply(metaDatas, names))
  metaDatas = lapply(metaDatas, function(metaData) {
    if ( all(columns %in% names(metaData)) ) return(metaData)
    missingColumns = lapply(setdiff(columns, names(metaData)), function(a) rep('', nrow(metaData)))
    names(missingColumns) = setdiff(columns, names(metaData))
    metaData = cbind(metaData, do.call(data.frame, missingColumns))
  })
  #merge and write to file
  metaData = do.call(rbind, metaDatas)
  write.table(metaData, file=targetMetaDataFile, sep='\t', quote=F, row.names=F)

  #merge data
  catLog('Loading data for merge')
  datas = mclapply(as.character(paths$Rdirectory), function(Rdirectory) {
    catLog('.')
    load(paste0(Rdirectory, '/allVariants.Rdata'))
    load(paste0(Rdirectory, '/clusters.Rdata'))
    allVariants = allVariants['variants']
    data = list(allVariants=allVariants, clusters=clusters)
    return(data)
  }, mc.cores=cpus)
  catLog('done.\n')
  #merge piece by piece
  variants = do.call(c, lapply(datas, function(d) d$allVariants$variants$variants))
  SNPs = do.call(rbind, lapply(datas, function(d) d$allVariants$variants$SNPs))
  SNPs = SNPs[!duplicated(rownames(SNPs)),]
  allVariants = list('variants'=list('variants'=variants, 'SNPs'=SNPs))
  ensureDirectoryExists(targetRdirectory)
  save(allVariants, file=variantsSaveFile)

  clusters = do.call(c, lapply(datas, function(d) d$clusters))
  save(clusters, file=clustersSaveFile)
}

#' looks for reccuring events across all individuals in a batch run
#'
#' @param metaDataFile character: path to the tab separated meta data file.
#' @param Rdirectory character: path to the Rdirectory of the batch run.
#' @param plotDirectory character: path to the plotDirectory of the batch run.
#' @param cpus integer: maximum number of threads used. default 1.
#' @param genome character: The genome build. 'hg19', hg38' or 'mm10'. Default 'hg19'
#' @param resourceDirectory character: path to the superFreq resources. Default 'superFreqResources'.
#' @param ...: Remaing argument are passed to superFreq::cohortAnalyseBatch()
#'
#' @details This function merges the data from multiple batches, setting them up for a joint cohort analysis. Note that this is only for the purpose of a downstream cohort analysis. The output R directory cannot be run as a batch Rdirectory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' superFreq(metaDataFile=metaDataFile,
#'           Rdirectory=Rdirectory,
#'           plotDirectory=plotDirectory,
#'           cpus=cpus,
#'           genome=genome)
#'
#' #The resource directory defaults to the same in superFreq() and superCohort()
#' #so do not need to be specified.
#' superCohort(metaDataFile=metaDataFile,
#'           Rdirectory=Rdirectory,
#'           plotDirectory=plotDirectory,
#'           cpus=cpus,
#'           genome=genome)
#'
#' }
superCohort = function(metaDataFile, Rdirectory='R', plotDirectory='plots', cpus=1, genome='hg19', resourceDirectory='superFreqResources', ...) {
  metaDataDir = paste0(dirname(metaDataFile), '/splitMetaData')
  metaDataFiles = list.files(metaDataDir, pattern='*.tsv', full.names=T)
  cohortDir = normalizePath(paste0(plotDirectory, '/cohort'))
  superFreq:::ensureDirectoryExists(cohortDir)
  batchDir = normalizePath(paste0(cohortDir, '/data'))
  superFreq:::ensureDirectoryExists(batchDir)
  Rdirectories = list.dirs(Rdirectory, recursive=F)
    
  paths = data.frame('metaDataFile' = metaDataFiles, 'Rdirectory'= Rdirectories)

  targetMetaDataFile = paste0(batchDir, '/metaData.tsv')
  targetRdirectory = paste0(batchDir, '/R')
  
  mergedPaths = mergeBatches(paths, targetMetaDataFile, targetRdirectory)
  
  metaData = importSampleMetaData(targetMetaDataFile)

  if ( !('PROJECT' %in% names(metaData)) )
    metaData$PROJECT = 'myProject'
  if ( !('GROUP' %in% names(metaData)) )
    metaData$GROUP = 'myGroup'
  
  cohortMetaDataFile = paste0(batchDir, '/cohortMetaData.tsv')
  write.table(metaData, file=cohortMetaDataFile, row.names=F, col.names=T, sep='\t', quote=F)
  
  outputDirectories = list('Rdirectory'=targetRdirectory)
  cohortAnalyseBatch(cohortMetaDataFile,
                     outputDirectories,
                     cpus=cpus,
                     genome=genome,
                     onlyDNA=F,
                     cosmicDirectory=paste0(resourceDirectory, '/COSMIC'),
                     ...)

}



checkRelatedness = function(Rdirectory='R', plotDirectory='plots', cpus=1) {
  samples = list.dirs(Rdirectory, recursive=F, full.names=F)
  variantFiles = paste0(Rdirectory, '/', samples, '/allVariants.Rdata')
  variantFiles = variantFiles[file.exists(variantFiles)]
  qsList = mclapply(variantFiles, function(file) {
        load(file)
        return(allVariants$variants$variants)
    }, mc.cores=cpus, mc.preschedule=F)
    qs = do.call(c, qsList)
    
    samples = names(qs)
    relatednessList = mclapply(samples, function(sample1) {
        cat(sample1, '\n')
        sapply(samples, function(sample2) {
            q1 = qs[[sample1]]
            q2 = qs[[sample2]]

            het1 = rownames(q1)[q1$db & !is.na(q1$dbMAF) & q1$dbMAF > 0.01 & q1$flag == '' &
              q1$var < 0.9*q1$cov & q1$var > 0.2*q1$cov & q1$cov > 20]
            het1loose = rownames(q1)[q1$var < 0.95*q1$cov & q1$var > 0.1*q1$cov & q1$cov > 5]
            het2 = rownames(q2)[q2$db & !is.na(q2$dbMAF) &  q2$dbMAF > 0.01 & q2$flag == '' &
              q2$var < 0.9*q2$cov & q2$var > 0.2*q2$cov & q2$cov > 20]
            het2loose = rownames(q2)[q2$var < 0.95*q2$cov & q2$var > 0.1*q2$cov & q2$cov > 5]

            agreeHet = (sum(het1 %in% het2loose) + sum(het2 %in% het1loose))/2
            hetScore = agreeHet/((length(het1) + length(het2))/2)

            return(hetScore)
        })
    }, mc.cores=5, mc.preschedule=F)
    names(relatednessList) = samples
    relatednessMx = do.call(rbind, relatednessList)
  
  cohortPlots = paste0(plotDirectory, '/cohortWide')
  superFreq:::ensureDirectoryExists(cohortPlots)
  pdf(paste0(cohortPlots, '/relatedness.pdf'), height=10, width=15)
  makeHeatmap(relatednessMx, Rowv=NA, Colv=NA, col='sunset')
  dev.off()

  return(invisible(relatednessMx))
}
