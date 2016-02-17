#########################################################################
#
# Functions to handle meta data
#
#########################################################################

isSaneMetaData = function(metaData) {
  #has all the right classes everywhere
  if ( class(metaData) != 'list' ) return('metaData not a list!')

  if ( !('samples' %in% names(metaData)) ) return('No entry in metadata called \'samples\'!')
  if ( !('normalSamples' %in% names(metaData)) ) return('No entry in metadata called \'normalSamples\'!')
  if ( !('individuals' %in% names(metaData)) ) return('No entry in metadata called \'individuals\'!')
  if ( !('project' %in% names(metaData)) ) return('No entry in metadata called \'project\'!')
  if ( !('paths' %in% names(metaData)) ) return('No entry in metadata called \'paths\'!')

  #information about each sample, and a separate list for the pool of normals
  if ( class(metaData$samples) != 'data.frame' ) return('metaData$samples is not a data.frame!')
  if ( class(metaData$normalSamples) != 'data.frame' ) return('metaData$normalSamples is not a data.frame!')

  sampleClass = c(
    'BAM'='character',
    'NAME'='character',
    'NORMAL'='logical',
    'INDIVIDUAL'='character',
    'CAPTUREREGIONS'='character',
    'PROJECT'='character',
    'BATCH'='character',
    'DATATYPE'='character',
    'RDIRECTORY'='character',
    'PLOTDIRECTORY'='character')

  for ( name in names(sampleClass) ) {
    if ( !(toupper(name) %in% names(metaData$samples)) )
      return(paste0('No column in metaData$samples called ', name ,'!'))
    if ( class(metaData$samples[[name]]) != sampleClass[name] )
      return(paste0('metaData$samples$', name,' is not of class ', sampleClass[name],'!'))

    if ( !(name %in% names(metaData$normalSamples)) )
      return(paste0('No column in metaData$normalSamples called ', name ,'!'))
    if ( class(metaData$normalSamples[[name]]) != sampleClass[name] )
      return(paste0('metaData$normalSamples$', name,' is not of class ', sampleClass[name],'!'))
  }

  #information about each individual, mainly where to save R data and where to plot.
  if ( class(metaData$individuals) != 'data.frame' ) return('metaData$individuals is not a data.frame!')
  sampleClass = c(
    'individual'='character',
    'Rdirectory'='character',
    'plotDirectory'='character')

  for ( name in names(sampleClass) ) {
    if ( !(name %in% names(metaData$individuals)) )
      return(paste0('No column in metaData$individuals called ', name ,'!'))
    if ( class(metaData$individuals[[name]]) != sampleClass[name] )
      return(paste0('metaData$individuals$', name,' is not of class ', sampleClass[name],'!'))
  }

  #information about each project, mainly where to save R data and where to plot.
  if ( class(metaData$project) != 'data.frame' ) return('metaData$project is not a data.frame!')
  sampleClass = c(
    'project'='character',
    'Rdirectory'='character',
    'plotDirectory'='character')

  for ( name in names(sampleClass) ) {
    if ( !(name %in% names(metaData$project)) )
      return(paste0('No column in metaData$project called ', name ,'!'))
    if ( class(metaData$project[[name]]) != sampleClass[name] )
      return(paste0('metaData$project$', name,' is not of class ', sampleClass[name],'!'))
  }

  #information about each project, mainly where to save R data and where to plot.
  if ( class(metaData$paths) != 'list' ) return('metaData$paths is not a list!')
  sampleClass = c(
    'sources'='character',
    'dataDirectory'='character')

  for ( name in names(sampleClass) ) {
    if ( !(name %in% names(metaData$paths)) )
      return(paste0('No column in metaData$paths called ', name ,'!'))
    if ( class(metaData$paths[[name]]) != sampleClass[name] )
      return(paste0('metaData$paths$', name,' is not of class ', sampleClass[name],'!'))
  }

  
  #paths point to existing files (where needed)
  bamFiles = c(metaData$samples$BAM, metaData$normalSamples$BAM)
  if ( any(!file.exists(bamFiles)) )
    return(paste0('The following bam files do not exist: ', do.call(paste, as.list(bamFiles[!file.exists(bamFiles)])), '.'))
  
  #unique row names

  #matching row and sample names

  #check that reference genomes exist for every sample

  return('yes')
}


makeMetaData = function(sourceFile, dataDirectory) {
  #load file
  samples = read.table(sourceFile, header=T, as.is=T, fill=T, sep=',')

  #set up all metadata entries, while sanity checking input
  metaData = list()

  if ( !('BAM' %in% names(samples)) ) stop('Must have a BAM column in the metadata$samples.')
  samples$BAM = gsub(' +$', '', samples$BAM)

  if ( !('NORMAL' %in% names(samples)) ) stop('Must have a NORMAL column in the metadata$samples.')
  samples$NORMAL = samples$NORMAL == 'YES'
  if ( !('NAME' %in% names(samples)) ) stop('Must have a NAME column in the metadata$samples.')
  samples$NAME = make.names(samples$NAME)
  if ( !('SAMPLE' %in% names(samples)) ) stop('Must have a SAMPLE column in the metadata$samples.')
  samples$SAMPLE = make.names(samples$SAMPLE)
  if ( ('PROJECT.SUBGROUP' %in% names(samples)) )
    samples$PROJECT.SUBGROUP[samples$PROJECT.SUBGROUP != ''] = make.names(samples$PROJECT.SUBGROUP)[samples$PROJECT.SUBGROUP != '']
  samples$RDIRECTORY = paste0(dataDirectory, '/R/samples/',samples$NAME)
  samples$PLOTDIRECTORY = paste0(dataDirectory, '/plots/samples/',samples$NAME)
  rownames(samples) = samples$NAME

  normalSamples = samples[samples$NORMAL,]
  metaData$normalSamples = normalSamples
  
  if ( !('INDIVIDUAL' %in% names(samples)) ) stop('Must have a INDIVIDUAL column in the metadata.')
  samples$INDIVIDUAL = make.names(samples$INDIVIDUAL)
  individuals = unique(samples$INDIVIDUAL)
  individuals = data.frame(
    'individual'=individuals,
    'Rdirectory'=paste0(dataDirectory, '/R/individuals/',individuals),
    'plotDirectory'=paste0(dataDirectory, '/plots/individuals/',individuals), stringsAsFactors=F)
  rownames(individuals) = individuals$individual
  metaData$individuals = individuals

  if ( !('PROJECT' %in% names(samples)) ) stop('Must have a PROJECT column in the metadata.')
  samples$PROJECT = make.names(samples$PROJECT)
  projects = unique(unlist(strsplit(samples$PROJECT, split='\\.')))
  projects = data.frame(
    'project'=projects,
    'Rdirectory'=paste0(dataDirectory, '/R/projects/',projects),
    'plotDirectory'=paste0(dataDirectory, '/plots/projects/',projects), stringsAsFactors=F)
  rownames(projects) = projects$project
  metaData$project = projects
  metaData$samples = samples


  paths = list(
    'sources'=sourceFile,
    'dataDirectory'=dataDirectory)
  metaData$paths = paths

  diagnostics = list(
    'samples'=lapply(metaData$samples$NAME, function(sample) list()),
    'individual'=lapply(metaData$individuals$individual, function(ind) list()),
    'project'=lapply(metaData$individuals$project, function(project) list()),
    'batch'=lapply(unique(metaData$samples$BATCH), function(batch) list()))
  names(diagnostics$samples) = metaData$samples$NAME
  names(diagnostics$individual) = metaData$individuals$individual
  names(diagnostics$project) = metaData$individuals$project
  names(diagnostics$batch) = unique(metaData$samples$BATCH)
  metaData$paths = paths

  #create directories if needed

  isSane = isSaneMetaData(metaData)
  if ( isSane == 'yes' )
    return(metaData)
  
  cat('Meta data not sane.\n', isSane, '\n')
  return(isSane)
}

createDirectories = function(metaData, linkCommand = function(sor, link) {paste0('ln -s ', sor, ' ', link)}) {
  isSaneMetaData(metaData)

  dataDirectory = metaData$paths$dataDirectory
  
  #level 1: data directory
  if ( !file.exists(metaData$paths$dataDirectory) ) stop('Data directory not found.')
  if ( file.exists(metaData$paths$dataDirectory) & !file.info(metaData$paths$dataDirectory)[,'isdir'] ) stop('Data directory already exists, but is not a directory. Please rename file or change data directory name.')

  #level 2: R and plots
  resourceDir = paste0(dataDirectory, '/resources')
  if ( !file.exists(resourceDir) ) dir.create(resourceDir)
  plotDir = paste0(dataDirectory, '/plots')
  if ( !file.exists(plotDir) ) dir.create(plotDir)
  RDir = paste0(dataDirectory, '/R')
  if ( !file.exists(RDir) ) dir.create(RDir)
  vcfDir = paste0(dataDirectory, '/vcf')
  if ( !file.exists(vcfDir) ) dir.create(vcfDir)
  bamDir = paste0(dataDirectory, '/bamLinks')
  if ( !file.exists(bamDir) ) dir.create(bamDir)

  #level 3: sample, individual and project subfolders
  dir = paste0(dataDirectory, '/plots/projects')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/plots/samples')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/plots/individuals')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/plots/normals')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/plots/byType')
  if ( !file.exists(dir) ) dir.create(dir)

  dir = paste0(dataDirectory, '/R/projects')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/R/samples')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/R/individuals')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/R/normals')
  if ( !file.exists(dir) ) dir.create(dir)

  dir = paste0(dataDirectory, '/bamLinks/projects')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/bamLinks/samples')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/bamLinks/individuals')
  if ( !file.exists(dir) ) dir.create(dir)
  dir = paste0(dataDirectory, '/bamLinks/normals')
  if ( !file.exists(dir) ) dir.create(dir)

  #level 4: create plot folders for each sample, individual, and project
  for ( project in getProjects(metaData, onlyDNA=F) ) {
    projectPlotDir = paste0(dataDirectory, '/plots/projects/', project)
    ensureDirectoryExists(projectPlotDir)
    projectRdir = paste0(dataDirectory, '/R/projects/', project)
    ensureDirectoryExists(projectRdir)
  }
  for ( individual in metaData$individuals$individual ) {
    individualPlotDir = paste0(dataDirectory, '/plots/individuals/', individual)
    ensureDirectoryExists(individualPlotDir)
    individualRdir = paste0(dataDirectory, '/R/individuals/', individual)
    ensureDirectoryExists(individualRdir)
  }
  for ( sample in metaData$samples$NAME ) {
    samplePlotDir = paste0(dataDirectory, '/plots/samples/', sample)
    ensureDirectoryExists(samplePlotDir)
    sampleRdir = paste0(dataDirectory, '/R/samples/', sample)
    ensureDirectoryExists(sampleRdir)
  }

  #level 5: link plot folders to subsets: projects to individuals, individuals to samples.
  for ( project in getProjects(metaData, onlyDNA=F) ) {
    individuals = unique(metaData$samples[inProject(metaData, project),]$INDIVIDUAL)
    for ( individual in individuals ) {
      sourceDir = paste0(dataDirectory, '/plots/individuals/', individual)
      linkDir = paste0(dataDirectory, '/plots/projects/', project, '/', individual)
      if ( file.exists(sourceDir) & !file.exists(linkDir) ) {
        cat(linkCommand(sourceDir, linkDir), '\n')
        system(linkCommand(sourceDir, linkDir), intern=T)
      }
    }
  }
  for ( individual in metaData$individuals$individual ) {
    samples = metaData$samples$NAME[metaData$samples$INDIVIDUAL == individual]
    for ( sample in samples ) {
      sourceDir = paste0(dataDirectory, '/plots/samples/', sample)
      linkDir = paste0(dataDirectory, '/plots/individuals/', individual, '/', sample)
      if ( file.exists(sourceDir) & !file.exists(linkDir) ) {
        cat(linkCommand(sourceDir, linkDir), '\n')
        system(linkCommand(sourceDir, linkDir), intern=T)
      }
    }
  }

  invisible('yes')
}

saveMetaData = function(metaData) {
  #save metadata to Rdata file, and print (relevant) information to csv

  #check if csv is already open etc, send warning if can't save.

  #fix permission
}


linkBams = function(metaData, linkCommand = function(sor, link) {paste0('ln -s ', sor, ' ', link)}) {
  isSaneMetaData(metaData)

  #samples
  for ( row in 1:nrow(metaData$samples) ) {
    sor = metaData$samples$BAM[row]
    link = paste0(metaData$paths$dataDirectory, '/bamLinks/samples/', metaData$samples$NAME[row], '.bam')
    if ( file.exists(sor) & !file.exists(link) ) {
      sor = gsub(' ', '\\\\ ', sor)
      link = gsub(' ', '\\\\ ', link)
      cat(linkCommand(sor, link), '\n')
      system(linkCommand(sor, link), intern=T)
    }
    else if ( !file.exists(sor) )
      cat(paste0(sor, ' does not exists!\n'))

    sor = paste0(metaData$samples$BAM[row], '.bai')
    link = paste0(metaData$paths$dataDirectory, '/bamLinks/samples/', metaData$samples$NAME[row], '.bam.bai')
    if ( file.exists(sor) & !file.exists(link) ) {
      sor = gsub(' ', '\\\\ ', sor)
      link = gsub(' ', '\\\\ ', link)
      cat(linkCommand(sor, link), '\n')
      system(linkCommand(sor, link), intern=T)
    }
    else if ( !file.exists(sor) )
      cat(paste0(sor, ' does not exists!\n'))
  }

  #individuals
  for ( ind in unique(metaData$samples$INDIVIDUAL) ) {
    is = which(metaData$samples$INDIVIDUAL == ind)
    indDir = paste0(metaData$paths$dataDirectory, '/bamLinks/individuals/', ind)
    if ( !file.exists(indDir) )
      dir.create(indDir)
    
    for ( row in is ) {
      sor = metaData$samples$BAM[row]
      link = paste0(indDir, '/', metaData$samples$NAME[row], '.bam')
      if ( file.exists(sor) & !file.exists(link) ) {
        sor = gsub(' ', '\\\\ ', sor)
        link = gsub(' ', '\\\\ ', link)
        cat(linkCommand(sor, link), '\n')
        system(linkCommand(sor, link), intern=T)
      }
      else if ( !file.exists(sor) )
        cat(paste0(sor, ' does not exists!\n'))
      
      sor = paste0(metaData$samples$BAM[row], '.bai')
      link = paste0(indDir, '/', metaData$samples$NAME[row], '.bam.bai')
      if ( file.exists(sor) & !file.exists(link) ) {
        sor = gsub(' ', '\\\\ ', sor)
        link = gsub(' ', '\\\\ ', link)
        cat(linkCommand(sor, link), '\n')
        system(linkCommand(sor, link), intern=T)
      }
      else if ( !file.exists(sor) )
        cat(paste0(sor, ' does not exists!\n'))
    }
  }

  #projects
  for ( project in unique(metaData$samples$PROJECT) ) {
    is = which(metaData$samples$PROJECT == project)
    proDir = paste0(metaData$paths$dataDirectory, '/bamLinks/projects/', project)
    if ( !file.exists(proDir) )
      dir.create(proDir)
    
    for ( row in is ) {
      sor = metaData$samples$BAM[row]
      link = paste0(proDir, '/', metaData$samples$NAME[row], '.bam')
      if ( file.exists(sor) & !file.exists(link) ) {
        sor = gsub(' ', '\\\\ ', sor)
        link = gsub(' ', '\\\\ ', link)
        cat(linkCommand(sor, link), '\n')
        system(linkCommand(sor, link), intern=T)
      }
      else if ( !file.exists(sor) )
        cat(paste0(sor, ' does not exists!\n'))
      
      sor = paste0(metaData$samples$BAM[row], '.bai')
      link = paste0(proDir, '/', metaData$samples$NAME[row], '.bam.bai')
      if ( file.exists(sor) & !file.exists(link) ) {
        sor = gsub(' ', '\\\\ ', sor)
        link = gsub(' ', '\\\\ ', link)
        cat(linkCommand(sor, link), '\n')
        system(linkCommand(sor, link), intern=T)
      }
      else if ( !file.exists(sor) )
        cat(paste0(sor, ' does not exists!\n'))
    }
  }

}

ensureDirectoryExists = function(file) {
  dir = dirname(file)
  if ( !file.exists(dir) ) dir.create(dir)
}


batchToAnalysisDirectory = function(batch) {
  warning('Called a hardcoding function: batchToAnalysisDirectory')
  if ( length(batch) > 1 ) sapply(batch, batchToAnalysisDirectory)

  if ( batch == 'P5_AGRF_CAGRF8961' ) return('/wehisan/general/academic/grp_leukemia_genomics/P5_AGRF_CAGRF8961_C5GJLACXX/analysis/')
  if ( batch == 'P6_AGRF_CAGRF9332' ) return('/wehisan/general/academic/grp_leukemia_genomics/P6_AGRF_CAGRF9332_C57VCANXX/')
  if ( batch == 'P7_AGRF_CAGRF9333' ) return('/wehisan/general/academic/grp_leukemia_genomics/P7_AGRF_CAGRF9333_C5M1HANXX/')
  if ( batch == 'P8plus_CAGRF9783' ) return('/wehisan/general/academic/grp_leukemia_genomics/P8plus_CAGRF9783/')
  if ( batch == 'P9_NextSeq_WEHI_20141223' ) return('/wehisan/general/academic/grp_leukemia_genomics/P9_NextSeq_WEHI_20141223/')
  if ( batch == 'P11plus_AGRF_CAGRF10495_C71R0ANXX' ) return('/wehisan/general/academic/grp_leukemia_genomics/P11plus_AGRF_CAGRF10495_C71R0ANXX/')
  if ( batch == 'F13TSFAPHT0461' ) return('/wehisan/general/academic/grp_leukemia_genomics/F13TSFAPHT0461_HUMwkbX/analysis-new/')
  if ( batch == 'AGRF_CAGRF8634' ) return('/wehisan/general/academic/grp_leukemia_genomics/AGRF_CAGRF8634_C4UJFACXX/analysis/')
  if ( batch == 'F13TSFAPHT0564' ) return('/wehisan/general/academic/grp_leukemia_genomics/')
  if ( batch == 'P4_AGRF_CAGRF8840' ) return('/wehisan/general/academic/grp_leukemia_genomics/')
  if ( batch == 'P3_AGRF_CAGRF8744' ) return('/wehisan/general/academic/grp_leukemia_genomics/')
  stop(paste0('Dont know about batch ', batch))
  return('')
}

