
#' Imports metadata about the samples
#'
#' @details This function imports metadata from a tab-separates file. The file should have a column name row, and then one row for each bam file. The metadata must include the following columns, but in any order.
#'
#'          BAM  Relative path to the bamfile from the metadata file.
#'          NAME The name of the sample that will appear on plots and in output. Must be unique.
#'          INDIVIDUAL Unique identifier for the individual. Samples from the same individual will be compared to each other, mutations will be tracked over samples of the individual, and a normal sample from the same individual will be used as matched normal.
#'          NORMAL If the sample is normal (tumor burden at most 1%). 'YES' or 'NO'
#'
#' @export
#' @examples
#' #import and print metadata
#' metaDataFile = '/absolute/path/to/myAnalysis/metaData.txt'
#' metaData = importMetaData(metaDataFile)
#' metaData
#'
importSampleMetaData = function(sampleMetaDataFile) {
  if ( !exists('catLog') ) assign('catLog', cat, envir=.GlobalEnv)
  if ( !file.exists(sampleMetaDataFile) ) stop("Meta data file ", sampleMetaDataFile, ' doesnt exist.')
  catLog('Loading sample meta data from file...')
  metaData = read.table(sampleMetaDataFile, header=T, as.is=T, fill=T, sep='\t', strip.white=T)
  if ( any(!(c('BAM', 'VCF', 'INDIVIDUAL', 'NAME', 'NORMAL', 'TIMEPOINT') %in% colnames(metaData))) )
    stop('Could not find required columns BAM, VCF, INDIVIDUAL, NAME, NORMAL, TIMEPOINT in sample meta data.\n
The meta data file should be a tab separated file with headings.\n')
  metaData$BAM = as.character(metaData$BAM)
  metaData$VCF = as.character(metaData$VCF)
  metaData$INDIVIDUAL = as.character(metaData$INDIVIDUAL)
  metaData$NAME = as.character(metaData$NAME)
  metaData$NORMAL = as.character(metaData$NORMAL)
  metaData$TIMEPOINT = as.character(metaData$TIMEPOINT)
  catLog('done.\n')
  newNames = make.names(c('normal', metaData$NAME), unique=T)[-1]
  if ( any(newNames != metaData$NAME) ) {
    catLog('Standardised sample names to:\n')
    catLog(newNames, sep='\n')
    catLog('\n')
  }
  metaData$NAME = newNames
  if ( any(duplicated(metaData$NAME)) ) {
    stop(paste('Duplicated NAME in meta data:', duplicated(metaData$NAME), '.\n'))
  }

  #resolve BAM path
  relativePath = !grepl('^[~/]', metaData$BAM)
  metaData$BAM = normalizePath(ifelse(relativePath,
    paste0(dirname(sampleMetaDataFile), '/', metaData$BAM),
    metaData$BAM))

  #resolve vcf path
  relativePath = !grepl('^[~/]', metaData$VCF)
  metaData$VCF = normalizePath(ifelse(relativePath,
    paste0(dirname(sampleMetaDataFile), '/', metaData$VCF),
    metaData$VCF))

  rownames(metaData) = metaData$NAME
  return(metaData)
}

#helper function to extract sample pairs from meta data
#cuts at 100 pairs.
metaToSamplePairs = function(names, individuals, normals) {
  catLog('Deciding which pairs to scatter plot..')
  pairs = list()
  for (individual in unique(individuals)) {
    rows = which(individual == individuals)
    if ( length(rows) < 2 ) next
    for ( row1 in rows ) {
      for ( row2 in rows[rows > row1] ) {
        if ( normals[row2] & !normals[row1] )
          pairs = c(pairs, list(c(names[row2], names[row1])))
        else
          pairs = c(pairs, list(c(names[row1], names[row2])))
      }
    }
  }
  if ( length(pairs) > 100 ) pairs = pairs[1:100]
  catLog('done.\n')
  return(pairs)
}

#helper function to extract time series from meta data.
metaToTimeSeries = function(names, individuals, normals) {
  catLog('Deciding which time series to plot..')
  series = list()
  for (individual in unique(individuals)) {
    rows = which(individual == individuals)
    if ( length(rows) < 1 ) next
    series = c(series, list(names[rows]))
    names(series)[length(series)] = individual
  }
  catLog('done.\n')
  return(series)
}
