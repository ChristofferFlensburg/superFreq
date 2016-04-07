#surely there must be a better way to track versions?
#' Return the current version of superFreq.
#' @details this is the version that goes in the logfile.
#'          First digit means very large changes.
#'          Second digit means algorithm changes that affects output.
#'          Third digit is debugging that shouldn't affect output of completed runs.
#'          1.0.0 will be the version used in the performance testing in the first preprint.
#' @export
superVersion = function() return('0.9.4')


#' Wrapper to run default superFreq analysis
#'
#' @param metaDataFile Character. A path to a tab separated files with headers:
#'                     BAM: path (absolute or relative from metaData file) to bam file of sample.
#'                          An index file is required as well.\cr
#'                     VCF: path (absolute or relative from metaData file) to vcf file of sample.
#'                          This should include both somatic and germline variants. The variants
#'                          will undergo further filtering, so it is not a problem if the .vcf
#'                          includes false positives, although large number of entries in the
#'                          .vcf increases run time.\cr
#'                     NAME: unique identifier of the sample. The names will be converted to
#'                           standard R names, which involves replacing separators (such as space,
#'                           dash or underscore) with dots.\cr
#'                     INDIVIDUAL: unique identifier of the indiviudal the sample comes from.
#'                                 This information is used pair up matched normal samples for
#'                                 somatic mutation filtering, and for identifying germline
#'                                 heterozygous SNPs. This column is also used to group samples
#'                                 that should be compared to each other: all samples from the
#'                                 same individual are analysed together, and mutations are
#'                                 tracked between the samples.\cr
#'                     NORMAL: Should be YES if the sample is normal, NO otherwise. A normal
#'                             sample is assumed to have no somatic mutations. In practice, a
#'                             tumor burden below a few percent works effectively as a normal,
#'                             but tumor burdens above 5% in a sample marked as normal can cause
#'                             severely reduced sensitivity in cancer samples from the same
#'                             indiviudal.\cr
#'                     TIMEPOINT: Mostly used as label in plots together with the sample. Can
#'                                be Diagnosis, Relapse, resistant or similar labels. Does not
#'                                influence the analyss otherwise. 
#' @param captureRegions Character. A path to a bed file of the unpadded capture regions.
#' @param normalDirectory Character. A path to a directory containing reference normals.
#'                        The directory should contain (links to) indexed .bam for the pool of
#'                        normal samples. These samples can, but dont have to, be matched normals
#'                        of the analysed samples. They need to be from the same capture as the
#'                        analysed samples, and preferably sequenced in the same lab. At least
#'                        two are required, but more are better (but slower). 10 is a good number.
#' @param Rdirectory Character. A path to the directory where data will be stored.
#'                   The directory will be created if needed, but the parent directory must exist.
#' @param plotDirectory Character. A path to the directory where plots and output is placed.
#'                      The directory will be created if needed, but the parent directory must exist.
#' @param reference Character. Path to the indexed fasta file the bams were aligned to.
#' @param genome Character. The reference genome the samples are aligned to, such as "hg19"
#'               or "mm10". Defaults to "hg19".
#' @param BQoffset Integer. The offset of the base quality encoding. Defaults to 33.
#' @param forceRedo named list of Booleans. Controls which step are retrieved from previously
#'                  saved data, and which are redone even if previous results are saved.
#'                  If a step is forced to be redone, depending downstream steps are also
#'                  redone. The named entries in the list should be:\cr
#'
#'                  forceRedoCount: The read counting over capture regions for the samples.\cr
#'                  forceRedoNormalCount: The read counting over capture regions for the
#'                                        set of reference normals.\cr
#'                  forceRedoFit: The differential coverage analysis.\cr
#'                  forceRedoVolcanoes: Volcano plots of the differential coverage.\cr
#'                  forceRedoDifferentRegions: Top table output of gained/lost genes/exons.\cr
#'                  forceRedoSNPs: Deprecated, usually not called.\cr
#'                  forceRedoVariants: Import and QC of the variants pointed to in the VCF. \cr
#'                  forceRedoNormalSNPs: Deprecated, usually not called.\cr
#'                  forceRedoNormalVariants: Import and QC of the variants in the set of
#'                                           reference normals.\cr
#'                  forceRedoMatchFlag: QC of variants using the set of reference normals.
#'                                      Somatic SNV are also separated from germline.\cr
#'                  forceRedoScatters: Scatter plots of the SNV frequencies between samples.\cr
#'                  forceRedoOutputSomatic: Top tables of somatic SNVs by sample.\cr
#'                  forceRedoNewVariants: Top table of new somatic SNVs compared to other
#'                                        samples.\cr
#'                  forceRedoSNPprogression: Heatmap of SNV frequencies over multiple samples.\cr
#'                  forceRedoCNV: Copy number calling.\cr
#'                  forceRedoCNVplots: Copy number plots.\cr
#'                  forceRedoSummary: Summary plot of SNVs and CNVs over all samples.\cr
#'                  forceRedoStories: Clonal tracking of somatic mutations.\cr
#'                  forceRedoRiver: Plots of the clonal evolution.\cr
#'                  forceRedoVEP: VEP call (Variant Effect Predictor) and COSMIC linking.
#'
#'                  The default value is forceRedoNothing(), setting all entries to FALSE.
#'                  To redo everything, set to forceRedoEverything().
#' @param cpus Integer. Maximum number of parallel threads. Often threaded over chromsomes,
#'             so numbers like 4, 6, 8, 12, 24 are good for human genomes. Default 3.
#' @param outputToTerminalAsWell Boolean. If the rather verbose log output should be sent to
#'                               terminal as well as to the logfile (in the Rdirectory).
#'                               Default TRUE.
#' @param systematicVariance Numeric. A measure of the expected systematic variance in coverage
#'                           log fold change. Larger value decreases false CNV calls, smaller
#'                           values increase sensitivity. Defaults to 0.03.
#' @param maxCov Integer. The coverage at which systematic variance of SNV frequencies are equal
#'                        to the Poissonian variance. Larger values increase sensitivity in CNV
#'                        calls, smaller values decrease false positives. Default 150.
#' @param dbSNPdirectory Character. The location of the directory where the superFreq dbSNP data
#'                       is located. If the directory doesn't exist, it will be created and the
#'                       data will be downloaded. Defaults to superFreqDbSNP (in the directory
#'                       where R is run).
#' @param cosmicDirectory Character. The location of the directory where the superFreq COSMIC data
#'                       is located. If the directory doesn't exist, it will be created and the
#'                       data will be downloaded. Defaults to superFreqCOSMIC (in the directory
#'                       where R is run).
#' @param mode Character. The mode to run in. The default 'Exome' is almost always used. An experimental
#                         'genome' is also allowed, that doesnt do MA correction for counts.
#'
#' @details This function runs a full SNV, SNP, CNV and clonality analysis in the input exome data.
#'          Note that a lot of settings go by default, such as base quality offset (33) and genome (hg19).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #minimal example to get it running, assuming data fits the defaults.
#' #you really want to skim through the other settings before running though.
#' metaDataFile = 'metaData.tsv'
#' captureRegions='captureRegions.bed'
#' normalDirectory = '../referenceNormals'
#' reference = '../reference/hg19.fa'
#' Rdirectory = 'R'
#' plotDirectory='plots'
#' 
#' data = superFreq(metaDataFile, captureRegions, normalDirectory,
#'                  Rdirectory, plotDirectory, reference)
#' }
superFreq = function(metaDataFile, captureRegions, normalDirectory, Rdirectory, plotDirectory, reference,
  genome='hg19', BQoffset=33, cpus=3, outputToTerminalAsWell=T, forceRedo=forceRedoNothing(),
  systematicVariance=0.03, maxCov=150, cloneDistanceCut=-qnorm(0.01), dbSNPdirectory='superFreqDbSNP', cosmicDirectory='superFreqCOSMIC', mode='exome') {

  inputFiles =
    superInputFiles(metaDataFile, captureRegions, normalDirectory, dbSNPdirectory,
                    reference)
  #sanityCheckSuperInputFiles(inputFiles)

  outputDirectories = superOutputDirectories(Rdirectory=Rdirectory, plotDirectory=plotDirectory)
  #sanityCheckSuperOutputDirectories(outputDirectories)

  settings = defaultSuperSettings(genome=genome, BQoffset=BQoffset)
  if ( mode == 'genome' ) settings$MAcorrection = F

  runtimeSettings = defaultSuperRuntimeSettings(cpus=cpus, outputToTerminalAsWell=outputToTerminalAsWell)

  parameters = defaultSuperParameters(systematicVariance=systematicVariance, maxCov=maxCov, cloneDistanceCut=cloneDistanceCut)

  downloadSuperFreqDbSNP(dbSNPdirectory)
  downloadSuperFreqCOSMIC(cosmicDirectory)

  
  analyse(inputFiles=inputFiles, outputDirectories=outputDirectories, settings=settings, forceRedo=forceRedo,
                 runtimeSettings=runtimeSettings, parameters=parameters, byIndividual=T)

  postAnalyseVEP(outputDirectories, inputFiles=inputFiles, cosmicDirectory=cosmicDirectory,
                 cpus=cpus, forceRedo=forceRedo$forceRedoVEP)

  
}





#' Analyse exomes
#'
#' @param inputFiles A named list of input files, containing the entries metaDataFile, vcfFiles, normalDirectory, captureRegionsFile and dbSNPdirectory
#' @param outputDirectories A named list of output directories, containing the entries Rdirectory and plotDirectory where the saved data and plots will be stored respectively.
#' @param settings A named list containing the entries genome and BQoffset. The only genome supporter atm is 'hg19', and the BQ offset is 33 for most exomes, altough some have 64. Check your fastqc files if you are not sure.
#' @param forceRedo A named list of logicals controling if existing saved data should be loaded or regenerated (overwriting the previous saved data). shortcuts to create these lists are forceRedoNothing() and forceRedoEverything().
#' @param runtimeSettings A named list containing the entries cpus and outputToTerminalAsWell. cpus is an integer controling the maximum number of cpus used in parallel, and outputToTerminalAsWell prints log data to the R session as well as to the log file.
#' @keywords analyse exomes CNV clonality
#'
#' @details This function runs a full SNV, SNP, CNV and clonality analysis in the input exome data.
#'          Set up input data as is shown in the example.
#'
#' @export
#'
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#'
#' @examples
#' \dontrun{
#' metaDataFile = '/absolute/path/to/metaData.txt'
#'
#' captureRegionsFile = '/absolute/path/to/captureRegions.bed'
#'
#' dbSNPdirectory = '/absolute/path/to/dbSNP'
#'
#' normalDirectory = '/path/to/normalDirectory'
#'
#' Rdirectory = '/absolute/path/to/R'
#' plotDirectory = '/absolute/path/to/plots'
#'
#' cpus=6
#' 
#' BQoffset = 33
#' genome = 'hg19'
#'
#'
#' inputFiles =
#'   list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFiles,
#'        'normalDirectory'=normalDirectory, 'normalCoverageDirectory'=normalCoverageDirectory,
#'        'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
#' outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
#' 
#' runtimeSettings = list('cpus'=cpus, 'outputToTerminalAsWell'=T)
#' settings = list('genome'=genome, 'BQoffset'=BQoffset)
#'
#' forceRedo = forceRedoNothing()
#'
#' parameters = list('systematicVariance'=0.03, 'maxCov'=150)
#'
#' data = analyse(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings, parameters)
#'
#' }
analyse = function(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings,
  parameters=defaultSuperParameters(), byIndividual=T) {
  #loadMethods(byIndividual=byIndividual)
  options(stringsAsFactors = F)
  options(scipen = 10)

  if ( !all(c('Rdirectory', 'plotDirectory') %in% names(outputDirectories)) )
    stop('outputDirectories need all entries: Rdirectory, plotDirectory.')
  Rdirectory = normalizePath(outputDirectories$Rdirectory)
  if ( !exists('Rdirectory') ) stop('Need to set Rdirectory.')
  if ( class(Rdirectory) != 'character' ) stop('Rdirectory needs to be of class character.')
  
  if ( !file.exists(Rdirectory) ) {
    dirSuccess = dir.create(Rdirectory)
    if ( !dirSuccess ) stop('Failed to create the Rdirectory ', Rdirectory,
                            '. Note that parent directory must exist.')
  }

  #set up logfile and log start of run.
  outputToTerminalAsWell = runtimeSettings$outputToTerminalAsWell
  if ( is.null(outputToTerminalAsWell) | class(outputToTerminalAsWell) != 'logical' ) outputToTerminalAsWell=F
  logFile = normalizePath(paste0(Rdirectory, '/runtimeTracking.log'))
  assign('catLog', function(...) cat(..., file=logFile, append=T), envir = .GlobalEnv)
  if ( outputToTerminalAsWell )
    assign('catLog', function(...) {cat(..., file=logFile, append=T); cat(...)}, envir = .GlobalEnv)

  catLog('\n\n\n', as.character(Sys.time()),
         '\n######################################################################\n')
  catLog('Running superFreq version', superVersion(), '\n')


  neededInput = c('metaDataFile', 'vcfFiles', 'normalDirectory', 'captureRegionsFile', 'dbSNPdirectory')
  if ( byIndividual ) neededInput = c('metaDataFile', 'normalDirectory', 'captureRegionsFile', 'dbSNPdirectory')
  if ( !all(neededInput %in% names(inputFiles)) )
    stop(paste(c('inputFiles need all entries:', neededInput), collapse=' '))
  if ( !all(c('BQoffset', 'genome') %in% names(settings)) )
    stop('settings need all entries: BQoffset, genome.')
  if ( !all(c('cpus') %in% names(runtimeSettings)) )
    stop('runtimeSettings need all entries: cpus.')
  
  cpus = runtimeSettings$cpus
  
  genome = settings$genome
  sampleMetaDataFile = inputFiles$metaDataFile
  vcfFiles = inputFiles$vcfFiles
  if ( 'reference' %in% names(inputFiles) ) reference = inputFiles$reference
  else reference = ''
  normalDirectory = inputFiles$normalDirectory
  normalRdirectory = paste0(normalDirectory, '/R')
  if ( !file.exists(normalRdirectory) ) {
    dirSuccess = dir.create(normalRdirectory)
    if ( !dirSuccess ) stop('Failed to create the normal R directory ', normalRdirectory,
                            '. Note that parent directory must exist.')
  }
  if ( 'normalCoverageDirectory' %in% names(inputFiles) )
    normalCoverageDirectory = inputFiles$normalCoverageDirectory
  else
    normalCoverageDirectory = inputFiles$normalDirectory
  normalCoverageRdirectory = paste0(normalCoverageDirectory, '/R')
  if ( !file.exists(normalCoverageRdirectory) ) {
    dirSuccess = dir.create(normalCoverageRdirectory)
    if ( !dirSuccess ) stop('Failed to create the normal coverage R directory ', normalCoverageRdirectory,
                            '. Note that parent directory must exist.')
  }
  dbSNPdirectory = inputFiles$dbSNPdirectory
  captureRegionsFile = inputFiles$captureRegionsFile

  plotDirectory = normalizePath(outputDirectories$plotDirectory)
  BQoffset = settings$BQoffset
  
  #sanity check input to the R script.
  if ( !exists('sampleMetaDataFile') ) stop('Need to set sampleMetaDataFile.')
  if ( class(sampleMetaDataFile) != 'character' ) stop('sampleMetaDataFile needs to be of class character.')
  if ( !file.exists(sampleMetaDataFile) ) stop('sampleMetaDataFile ', sampleMetaDataFile, ' not found.')

  if ( !byIndividual ) {
    if ( !exists('vcfFiles') ) stop('Need to set vcfFiles!')
    if ( class(vcfFiles) != 'character' ) stop('vcfFiles needs to be of class character.')
    for ( vcfFile in vcfFiles ) if ( !file.exists(vcfFile) ) stop('vcfFile ', vcfFile, ' not found.')
  }
  
  if ( !exists('normalDirectory') ) stop('Need to set normalDirectory!')
  if ( class(normalDirectory) != 'character' ) stop('normalDirectory needs to be of class character.')
  if ( !file.exists(normalDirectory) ) stop('normalDirectory ', normalDirectory, ' not found.')

    if ( !exists('normalCoverageDirectory') ) stop('Need to set normalCoverageDirectory!')
  if ( class(normalCoverageDirectory) != 'character' ) stop('normalCoverageDirectory needs to be of class character.')
  if ( !file.exists(normalCoverageDirectory) ) stop('normalCoverageDirectory ', normalCoverageDirectory, ' not found.')
    
  if ( !exists('plotDirectory') ) stop('Need to set plotDirectory.')
  if ( class(plotDirectory) != 'character' ) stop('plotDirectory needs to be of class character.')
  if ( !file.exists(plotDirectory) ) {
    dirSuccess = dir.create(plotDirectory)
    if ( !dirSuccess ) stop('Failed to create the plotDirectory ', plotDirectory,
                            '. Note that parent directory must exist.')
  }



  cat('Runtime tracking and QC information printed to ', logFile, '.\n', sep='')
  catLog('Starting run with input files:',
         '\nsampleMetaDataFile:', sampleMetaDataFile,
         '\nvcfFiles:\n')
  catLog(vcfFiles, sep='\n')
  catLog('Normal directory:', normalDirectory, '\n')
  catLog('Normal coverage directory:', normalCoverageDirectory, '\n')
  catLog('dbSNP directory:', dbSNPdirectory, '\n')
  catLog('capture regions:', captureRegionsFile, '\n')
  catLog('Plotting to', plotDirectory, '\n')
  catLog('Saving R files to', Rdirectory, '\n')
  catLog('Genome is', genome, '\n')
  catLog('Running on at most', cpus, 'cpus.\n')


  if ( !(genome %in% c('hg19', 'mm10')) ) stop('Only genomes that are supported atm are hg19 and mm10, sorry.\nNew genomes can easily be added though, please contact the authors.\n')

  if ( class(parameters) != 'list' ) stop('parameters need to be of class list. Dont provide this to analyse() for default settings, otherwise a named list with the parameters.')
  catLog('\nParameters for this run are:\n')
  if ( 'maxCov' %in% names(parameters) &  class(parameters$maxCov) != 'numeric' ) stop('parameter maxCov needs to be numeric.')
  if ( 'maxCov' %in% names(parameters) )   assign('.maxCov', parameters$maxCov, envir = .GlobalEnv)
  else assign('.maxCov', 150, envir = .GlobalEnv)
  catLog('   maxCov:              ', get('.maxCov', envir = .GlobalEnv), '\n', sep='')
  if ( 'systematicVariance' %in% names(parameters) &  class(parameters$systematicVariance) != 'numeric' ) stop('parameter systematicVariance needs to be numeric.')
  if ( 'systematicVariance' %in% names(parameters) ) assign('.systematicVariance', parameters$systematicVariance, envir = .GlobalEnv)
  else assign('.systematicVariance', 0.03, envir = .GlobalEnv)
  catLog('   systematicVariance:  ', get('.systematicVariance', envir = .GlobalEnv), '\n', sep='')

  if ( 'cloneDistanceCut' %in% names(parameters) &  class(parameters$cloneDistanceCut) != 'numeric' ) stop('parameter cloneDistanceCut needs to be numeric.')
  if ( 'cloneDistanceCut' %in% names(parameters) ) assign('.cloneDistanceCut', parameters$cloneDistanceCut, envir = .GlobalEnv)
  else assign('.cloneDistanceCut', -qnorm(0.01), envir = .GlobalEnv)
  catLog('   cloneDistanceCut:  ', get('.cloneDistanceCut', envir = .GlobalEnv), '\n', sep='')


  catLog('\n')

  #set forceRedo parameters to false unless already specified
  forceRedoCount = forceRedo$forceRedoCount
  forceRedoNormalCount = forceRedo$forceRedoNormalCount
  forceRedoFit = forceRedo$forceRedoFit
  forceRedoVolcanoes = forceRedo$forceRedoVolcanoes
  forceRedoDifferentRegions = forceRedo$forceRedoDifferentRegions
  forceRedoSNPs = forceRedo$forceRedoSNPs
  forceRedoVariants = forceRedo$forceRedoVariants
  forceRedoNormalSNPs = forceRedo$forceRedoNormalSNPs
  forceRedoNormalVariants = forceRedo$forceRedoNormalVariants
  forceRedoMatchFlag = forceRedo$forceRedoMatchFlag
  forceRedoScatters = forceRedo$forceRedoScatters
  forceRedoOutputSomatic = forceRedo$forceRedoOutputSomatic
  forceRedoNewVariants = forceRedo$forceRedoNewVariants
  forceRedoSNPprogression = forceRedo$forceRedoSNPprogression
  forceRedoCNV = forceRedo$forceRedoCNV
  forceRedoCNVplots = forceRedo$forceRedoCNVplots
  forceRedoSummary = forceRedo$forceRedoSummary
  forceRedoStories = forceRedo$forceRedoStories
  forceRedoRiver = forceRedo$forceRedoRiver

#make sure that if something is redone, then all depending steps are redone as well
  if ( forceRedoCount ) {
    catLog('Redoing coverage counts, so need to redo linear analysis.\n')
    forceRedoFit = T
  }
  if ( forceRedoNormalCount ) {
    catLog('Redoing normal coverage counts, so need to redo linear analysis of coverage.\n')
    forceRedoFit = T
  }
  if ( forceRedoFit ) {
    catLog('Redoing linear analysis of coverage, so need to redo CNVs, volcano plots and diffent regions sheet.\n')
    forceRedoCNV = T
    forceRedoVolcanoes = T
    forceRedoDifferentRegions = T
  }
  if ( forceRedoSNPs ) {
    catLog('Redoing SNPs, so need to redo quality flagging of variants.\n')
    forceRedoVariants = T
  }
  if ( forceRedoVariants ) {
    catLog('Redoing quality flagging of variants, so need to redo flagmatching with normals.\n')
    forceRedoMatchFlag = T
  }
  if ( forceRedoNormalSNPs ) {
    catLog('Redoing normal SNPs, so need to redo quality flagging of normal variants.\n')
    forceRedoNormalVariants = T
  }
  if ( forceRedoNormalVariants ) {
    catLog('Redoing normal quality flagging of variants, so need to redo flag matching with normals.\n')
    forceRedoMatchFlag = T
  }
  if ( forceRedoMatchFlag ) {
    catLog('Redoing flagmatching of normals, so need to redo frequency scatters, variant sheets, frequency progressions and CNVs.\n')
    forceRedoSNPprogression = T
    forceRedoCNV = T
  }
  if ( forceRedoCNV ) {
    catLog('Redoing CNVs, so need to redo CNV plots, summary plot and clonality stories.\n')
    forceRedoCNVplots = T
    forceRedoSummary = T
    forceRedoStories = T
  }
  if ( forceRedoStories ) {
    catLog('Redoing stories, so need to redo river and scatter plots, and new and somatic spread sheets.\n')
    forceRedoRiver = T
    forceRedoScatters = T
    forceRedoNewVariants = T
    forceRedoOutputSomatic = T
  }

  externalNormalBams =
    c(list.files(path=paste0(normalDirectory), pattern = '*.bam$', full.names=T),
      list.files(path=paste0(normalDirectory, '/bam'), pattern = '*.bam$', full.names=T))
  names(externalNormalBams) = gsub('.bam$', '', basename(externalNormalBams))
  catLog('Normal bamfiles are:\n')
  catLog(externalNormalBams, sep='\n')

  externalNormalCoverageBams =
        c(list.files(path=paste0(normalCoverageDirectory), pattern = '*.bam$', full.names=T),
      list.files(path=paste0(normalCoverageDirectory, '/bam'), pattern = '*.bam$', full.names=T))
  names(externalNormalCoverageBams) = gsub('.bam$', '', basename(externalNormalCoverageBams))
  catLog('Normal coverage bamfiles are:\n')
  catLog(externalNormalCoverageBams, sep='\n')

  captureRegions = try(importCaptureRegions(captureRegionsFile, reference=reference, Rdirectory=Rdirectory, genome=genome))
  if ( class(captureRegions) != 'GRanges' ) {
    catLog('Failed to import capture regions, aborting.\n')
    stop('Failed to import capture regions, aborting.\n')
  }
  if ( length(captureRegions) == 0 ) {
    catLog('Empty capture regions, aborting.\n')
    stop('Empty capture regions, aborting.\n')
  }
  if ( !('region' %in% colnames(mcols(captureRegions)) & 'gc' %in% colnames(mcols(captureRegions))) ) {
    catLog('Need both region name and gc content of capture regions, aborting.\n')
    stop('Need both region name and gc content of capture regions, aborting.\n')
  }
  if ( any(is.na(captureRegions$gc)) ) {
    catLog('NA gc content value in capture regions, aborting.\n')
    stop('NA gc content value in capture regions, aborting.\n')
  }
  catLog('Imported capture regions with', length(captureRegions), 'regions and',
         length(unique(captureRegions$region)), 'unique gene names.\n')
  catLog('Mean GC content is ', round(mean(captureRegions$gc), 3), '.\n', sep='')
  
  #read in metadata
  sampleMetaData = try(importSampleMetaData(sampleMetaDataFile))
  if ( class(sampleMetaData) != 'data.frame' ) {
    catLog('Failed to import meta data, aborting.\n')
    stop('Failed to import meta data, aborting.\n')
  }
  if ( !all(c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL') %in% colnames(sampleMetaData)) ) {
    missing = c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL')[!(c('BAM', 'INDIVIDUAL', 'NAME', 'TIMEPOINT', 'NORMAL') %in% colnames(sampleMetaData))]
    catLog('Missing columns in meta data:' , missing, ', aborting.\n')
    stop('Missing columns in meta data:' , missing, ', aborting.\n')
  }
  if ( !all(sampleMetaData$NORMAL %in% c('YES', 'NO')) ) {
    catLog('Want only YES or NO in normal column, aborting.\n')
    stop('Want only YES or NO in normal column, aborting.\n')
  }
  bamFiles = sampleMetaData$BAM
  names = make.names(sampleMetaData$NAME, unique=T)
  individuals = sampleMetaData$INDIVIDUAL
  timePoints = sampleMetaData$TIMEPOINT
  names(timePoints) = names(individuals) = names
  normals = as.logical(gsub('YES', 'T', gsub('NO', 'F', sampleMetaData$NORMAL)))
  names(normals) = names
  samplePairs = metaToSamplePairs(names, individuals, normals)
  timeSeries = metaToTimeSeries(names, individuals, normals)

  if ( any(!file.exists(bamFiles)) ) {
    catLog('Missing (or misnamed) bam files:' , bamFiles[!file.exists(bamFiles)], '. This will crash the run unless coverage and variants are already imported in previously saved data.\n')
    warning('Missing (or misnamed) bam files:' , bamFiles[!file.exists(bamFiles)], '. This will crash the run unless coverage and variants are already imported in previously saved data.\n')
  }
  bamIndexFiles = paste0(bamFiles, '.bai')
  bamIndexFiles2 = gsub('.bam$', '.bai', bamFiles)
  if ( any(!(file.exists(bamIndexFiles) | file.exists(bamIndexFiles2))) ) {
    missingIndex = !(file.exists(bamIndexFiles) | file.exists(bamIndexFiles2))
    catLog('Could not find bam index files for:' , bamFiles[missingIndex], '. This will crash the run unless coverage and variants are already imported in previously saved data.\n')
    warning('Could not find bam index files for:' , bamFiles[missingIndex], '. This will crash the run unless coverage and variants are already imported in previously saved data.\n')
  }
  
  catLog('##################################################################################################\n\n',
         as.character(Sys.time()),'\n',
         'Imported and sanity checked meta data. Looking good so far!\n',
         'metadata:\n')
  catLog('',colnames(sampleMetaData), '\n', sep='   ')
  for ( row in 1:nrow(sampleMetaData) )
    catLog('',as.matrix(sampleMetaData[row,]), '\n', sep='   ')
  catLog('\n timeSeries:\n')
  for ( ts in timeSeries )
    catLog('',ts, '\n', sep='   ')
  catLog('\n##################################################################################################\n\n')

  #compare coverage of samples to the pool of normals, using limma-voom.
  fit = try(runDE(bamFiles, names, externalNormalCoverageBams, captureRegions, Rdirectory, plotDirectory, genome=genome,
    normalCoverageRdirectory, settings=settings, cpus=cpus, forceRedoFit=forceRedoFit, forceRedoCount=forceRedoCount,
    forceRedoNormalCount=forceRedoNormalCount))
  if ( class(fit) == 'try-error' ) {
    catLog('Error in runDEforSamples! Input was:')
    catLog('bamFiles:', bamFiles, '\nnames:', names, '\nexternalNormalCoverageBams:', externalNormalCoverageBams,
           '\nstart(captureRegions)[1:4]:', start(captureRegions)[1:4], '\nRdirectory:', Rdirectory,
           '\nplotDirectory:', plotDirectory, '\nnormalCoverageRdirectory:', normalCoverageRdirectory, '\ncpus:', cpus,
           '\nforceRedoFit:', forceRedoFit, '\nforceRedoCount', forceRedoCount,
           '\nforceRedoNormalCount', forceRedoNormalCount, '\n')
    dumpInput(Rdirectory, list('bamFiles'=bamFiles, 'names'=names, 'externalNormalCoverageBams'=externalNormalCoverageBams,
                               'captureRegions'=captureRegions, 'Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory,
                               'normalCoverageRdirectory'=normalCoverageRdirectory, 'cpus'=cpus, 'forceRedoFit'=forceRedoFit,
                               'forceRedoCount'=forceRedoCount, 'forceRedoNormalCount'=forceRedoNormalCount))
    stop('Error in runDEforSamples.')
  }

  #Plot volcanoes and output an excel file with top DE regions.
  ret = try(makeFitPlots(fit, plotDirectory, genome,
    forceRedoVolcanoes=forceRedoVolcanoes, forceRedoDifferentRegions=forceRedoDifferentRegions))
  if ( class(ret) == 'try-error' ) {
    catLog('Error in makeFitPlots, will continue analysis anyway.')
    warning('Error in makeFitPlots.')
  }


  saveFile = paste0(Rdirectory, '/allVariants.Rdata')
  {
    if ( file.exists(saveFile) & !forceRedoMatchFlag & !forceRedoVariants & !forceRedoNormalVariants ) {
      catLog('Loading final version of combined variants.\n')
      load(file=saveFile)
    }
    else {
      #import, filter and QC the variants. Save to file.
      #The information about normals is used for QC, as there will be only true frequencies of 0, 0.5 and 1 in those samples.
      if ( byIndividual )
        variants = try(getVariantsByIndividual(sampleMetaData, captureRegions, genome, BQoffset, dbSNPdirectory,
          Rdirectory, plotDirectory, cpus=cpus, forceRedo=forceRedoVariants))
      else
        variants = try(getVariants(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset, dbSNPdirectory,
          Rdirectory, plotDirectory, cpus=cpus, forceRedoSNPs=forceRedoSNPs,
          forceRedoVariants=forceRedoVariants))
      if ( class(variants) == 'try-error' ) {
        catLog('Error in getVariants.\n')
        stop('Error in getVariants.')
      }
      
      #Get variants from the external normals
      normalVariants =
        try(getNormalVariants(variants, externalNormalBams, names(externalNormalBams), captureRegions,
                        genome, BQoffset, dbSNPdirectory, normalRdirectory, Rdirectory, plotDirectory, cpus=cpus,
                        forceRedoSNPs=forceRedoNormalSNPs, forceRedoVariants=forceRedoNormalVariants))
      if ( class(normalVariants) == 'try-error' ) {
        catLog('Error in getVariants for normals.\n')
        dumpInput(Rdirectory, list('variants'=variants, 'externalNormalBams'=externalNormalBams,
                                   'names'=names(externalNormalBams), 'captureRegions'=captureRegions,
                                   'genome'=genome, 'BQoffset'=BQoffset, 'Rdirectory'=Rdirectory,
                                   'normalRdirectory'=normalRdirectory, 'plotDirectory'=plotDirectory, 'cpus'=cpus,
                                   'forceRedoSNPs'=forceRedoNormalSNPs, 'forceRedoVariants'=forceRedoNormalVariants))
        stop('Error in getVariants for normals.')
      }
      
      #share variants with normals
      flaggingVersion = 'old'
      if ( 'flaggingVersion' %in% names(settings) ) flaggingVersion = settings$flaggingVersion
      allVariants = try(matchFlagVariants(variants, normalVariants, individuals, normals, genome,
        Rdirectory, flaggingVersion=flaggingVersion, cpus=cpus, byIndividual=byIndividual,
        forceRedoMatchFlag=forceRedoMatchFlag))
      if ( class(allVariants) == 'try-error' | !all(c('variants', 'normalVariants') %in% names(allVariants)) ) {
        catLog('Error in matchFlagVariants.\n')
        dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'individuals'=individuals, 'normals'=normals,
                                   'Rdirectory'=Rdirectory, 'forceRedoMatchFlag'=forceRedoMatchFlag))
        stop('Error in matchFlagVariants.')
      }
    }
  }
  variants = allVariants$variants
  normalVariants = allVariants$normalVariants
  a=try(setVariantLoss(normalVariants$variants))
  if ( class(a) == 'try-error' ) {
    catLog('Error in setVariantLoss(normalVariants)\n')
    stop('Error in setVariantLoss(normalVariants)!')
  }
    

  #call CNVs compared to the normals.
  cnvs =
    try(callCNVs(variants=variants, normalVariants=normalVariants, fitS=fit$fit, variants$SNPs,
                 names=names, individuals=individuals, normals=normals, Rdirectory=Rdirectory, plotDirectory=plotDirectory,
                 genome=genome, cpus=cpus, forceRedoCNV=forceRedoCNV))
  if ( class(cnvs) == 'try-error' ) {
    catLog('Error in callCNVs!\n')
    dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'fitS'=fit$fit,
                               'names'=names, 'individuals'=individuals, 'normals'=normals, 'Rdirectory'=Rdirectory,
                               'genome'=genome, 'cpus'=cpus, 'forceRedoCNV'=forceRedoCNV))
    stop('Error in callCNVs.')
  }


  #make CNV plots
  cnvplot = try(makeCNVplots(cnvs, plotDirectory=plotDirectory, genome, forceRedoCNVplots=forceRedoCNVplots))
  if ( class(cnvplot) == 'try-error' ) {
    catLog('Error in makeCNVplots! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeCNVplots! Continuing, but these plots are kindof useful.')
  }

  #combine SNPs and CNVs into stories of subclones.
  stories = try(getStories(variants=variants, normalVariants=normalVariants, cnvs=cnvs, timeSeries=timeSeries, normals=normals, genome=genome, cloneDistanceCut=get('.cloneDistanceCut', envir = .GlobalEnv), Rdirectory=Rdirectory,
    plotDirectory=plotDirectory, cpus=cpus, forceRedo=forceRedoStories))
  if ( class(stories) == 'try-error' ) {
    catLog('Error in getStories!\n')
    dumpInput(Rdirectory, list('variants'=variants, 'normalVariants'=normalVariants, 'cnvs'=cnvs,
                               'timeSeries'=timeSeries, 'Rdirectory'=Rdirectory,
                               'plotDirectory'=plotDirectory, 'cpus'=cpus, 'forceRedo'=forceRedoStories))
    stop('Error in getStories!')
  }
  variants = stories$variants
  normalVariants = stories$normalVariants
  stories = stories$stories

  #do multi-sample heatmaps and frequency progression
  progression = try(makeSNPprogressionPlots(variants, timeSeries, normals, plotDirectory, cpus=cpus, forceRedo=forceRedoSNPprogression))
  if ( class(progression) == 'try-error' ) {
    catLog('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeSNPprogressionPlots! Continuing anyway, but these plots are kindof useful.')
  }
  
  #make summary plots
  summary = try(makeSummaryPlot(variants, cnvs, normals, individuals, timePoints, plotDirectory,
    genome, cpus, forceRedo=forceRedoSummary))
  if ( class(summary) == 'try-error' ) {
    catLog('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.\n')
    warning('Error in makeSummaryPlot! Continuing, but these plots are kindof useful.')
  }
  
  #identify and output new variants
  newVar = try(outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, forceRedo=forceRedoNewVariants))
  if ( class(newVar) == 'try-error' ) {
    catLog('Error in outputNewVariants! Continuing anyway.\n')
    warning('Error in outputNewVariants! Continuing anyway.')
  }
  
  #output somatic variants
  somatics = try(outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, forceRedo=forceRedoOutputSomatic))
  if ( class(somatics) == 'try-error' ) {
    catLog('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.\n')
    warning('Error in outputSomaticVariants! Continuing anyway, but this is fairly important information.')
  }

  river = try(makeRiverPlots(stories, variants, genome, cpus=cpus, plotDirectory, forceRedo=forceRedoRiver))
  if ( class(river) == 'try-error' ) {
    catLog('Error in makeRiverPlots!\n')
    warning('Error in makeRiverPlots!')
  }

  #Make the scatter plots of the pairs
  scatter = try(makeScatterPlots(variants, samplePairs, timePoints, plotDirectory,
    genome=genome, cpus=cpus, forceRedo=forceRedoScatters))
  if ( class(scatter) == 'try-error' ) {
    catLog('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.\n')
    warning('Error in makeScatterPlots! Continuing anyway, but these plots are kindof useful.')
  }

  catLog('Run done! Have fun with the output! :)\n\n')
  
  return(list('fit'=fit, 'variants'=variants, 'normalVariants'=normalVariants, 'cnvs'=cnvs, 'stories'=stories))
}

#' Loads saved data
#'
#' This function returns the data produced from an analyse() run.
#' @param Rdirectory A character string pointing to the Rdirectory of the analysis.
#' @param setVariantLoss If the average variant loss should be calculated. This measures reference bias of the variant frequencies, is set globally and is used by many other functions.
#' @keywords load saved data
#' @export
#' @examples
#' \dontrun{
#' not run
#' data = loadData('~/myAnalysis/R', setVariantLoss=T)
#' cnvs = data$clusters
#' plotCR(cnvs[['aSample']]$clusters)
#' }
loadData = function(Rdirectory, setVariantLoss=F) {
  saveFiles = list.files(Rdirectory, pattern = '*.Rdata$', full.names=T)
  names = gsub('.Rdata$', '', basename(saveFiles))
  #if allVariants is done, some of the earlier data doesnt have to be loaded.
  if ( 'allVariants' %in% names ) {
    dontLoad =
      names %in% c('variants', 'normalVariants', 'SNPs', 'allVariantsPreVEP', 'variantsBI', 'normalVariantsBI') |
    grepl('positions\\..+', names) | grepl('variants\\..+', names) | grepl('normalVariants\\..+', names)
    
    saveFiles = saveFiles[!dontLoad]
    names = names[!dontLoad]
  }
  if ( 'fit' %in% names ) {
    saveFiles = saveFiles[!(names %in% c('fitP'))]
    names = names[!(names %in% c('fitP'))]
  }
  names(names) = saveFiles
  cat('Loading..')
  for ( file in saveFiles ) {
    cat(gsub('.Rdata$', '', basename(file)), '..', sep='')
    names[file] = load(file=file)
  }
  cat('done.\n')
  ret = lapply(names, function(name) get(name))
  names(ret) = gsub('.Rdata$', '', basename(saveFiles))

  if ( setVariantLoss & 'allvariants' %in% names(ret) )
    setVariantLoss(ret$allVariants$normalVariants$variants)

  return(ret)
}

#' Not needed in the package, ignore this function.
#'
#' This function loads the analysis functions used in analyse().
#' @keywords load methods analyse
#'
#' @examples
#' loadMethods()
loadMethods = function(stringsAsFactors = FALSE, byIndividual=F) {
  options(stringsAsFactors = stringsAsFactors)
  options(scipen = 10)
  assign('catLog', function(...) cat(...), envir = .GlobalEnv)

  libraryLoaded = library(limma, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load limma.\n')
    stop('Failed to load limma.')
  }
  libraryLoaded = library(edgeR, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load edgeR.\n')
    stop('Failed to load edgeR.')
  }
  libraryLoaded = library(Rsubread, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load Rsubread.\n')
    stop('Failed to load Rsubread.')
  }
  libraryLoaded = library(parallel, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load parallel.\n')
    stop('Failed to load parallel.')
  }
  libraryLoaded = library(Rsamtools, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load Rsamtools.\n')
    stop('Failed to load Rsamtools.')
  }
  libraryLoaded = library(GenomicRanges, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load GenomicRanges.\n')
    stop('Failed to load GenomicRanges.')
  }
  libraryLoaded = library(R.oo, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load R.oo.\n')
    stop('Failed to load R.oo.')
  }
  libraryLoaded = library(rtracklayer, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load rtracklayer.\n')
    stop('Failed to load rtracklayer.')
  }
  libraryLoaded = library(WriteXLS, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load WriteXLS.\n')
    stop('Failed to load WriteXLS.')
  }
  libraryLoaded = library(fields, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load fields.\n')
    stop('Failed to load fields.')
  }
  libraryLoaded = library(seqinr, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load seqinr.\n')
    stop('Failed to load seqinr.')
  }
  libraryLoaded = library(biomaRt, logical.return=T)
  if ( !libraryLoaded ) {
    catLog('Error: Failed to load biomaRt.\n')
    stop('Failed to load biomaRt.')
  }
  catLog('All libraries loaded successfully.\n')
  
  source('debug.R')
  source('addBindingStrength.R')
  source('runVEP.R')
  source('importCaptureRegions.R')
  source('importSampleMetaData.R')
  source('runDEexon.R')
  source('XRank.R')
  source('makeFitPlots.R')
  if ( byIndividual ) {
    catLog('Running in by-individual mode.\n')
    source('getVariantsByIndividual.R')
  }
  else
    source('getVariants.R')
  source('matchFlagVariants.R')
  source('makeScatterPlots.R')
  source('outputNewVariants.R')
  source('outputSomaticVariants.R')
  source('makeSNPprogressionPlots.R')
  source('callCNVs.R')
  source('CNVnormalisation.R')
  source('makeCNVplots.R')
  source('makeSummaryPlot.R')
  source('getStories.R')
  source('makeRiverPlots.R')

  if ( !exists('.maxCov', envir = .GlobalEnv) ) {
    catLog('Setting maxCov to default 150.\n')
    assign('.maxCov', 150, envir = .GlobalEnv)
  }
  if ( !exists('.systematicVariance', envir = .GlobalEnv) ) {
    catLog('Setting systematicVariance to default 0.03.\n')
    assign('.systematicVariance', 0.03, envir = .GlobalEnv)
  }
}

#' returns input that uses saved data if present.
#'
#' @details This function returns the input 'forceRedo' for analyse(), so that saved data
#'          from previous runs is used if present. This should be the input normally.
#' @export
#' @examples
#' \dontrun{
#' cpus=4
#' metaDataFile = '/absolute/path/to/myAnalysis/metaData.txt'
#' vcfFiles = list.files('/absolute/path/to/myAnalysis/vcf', pattern='chr[0-9XYMT]+.vcf$', full.names=T)
#' captureRegionsFile = '/absolute/path/to/myAnalysis/captureRegions.gc.bed'
#' dbSNPdirectory = '/absolute/path/to/myAnalysis/dbSNP'
#' normalDirectory = '/absolute/path/to/myAnalysis/normal'
#' Rdirectory = '/absolute/path/to/myAnalysis/R/'
#' plotDirectory = '/absolute/path/to/myAnalysis/plots/'
#' 
#' #The base quality phred offset. This can be read from fastqc analysis for example.
#' BQoffset = 33
#' genome = 'hg19'
#' 
#' inputFiles = list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFiles, 'normalDirectory'=normalDirectory,
#' 'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
#' outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
#' runtimeSettings = list('cpus'=cpus, 'outputToTerminalAsWell'=T)
#'
#' forceRedo = forceRedoNothing()
#' parameters = list('systematicVariance'=0.03, 'maxCov'=150)
#'
#' data = analyse(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings, parameters=parameters)
#' }
forceRedoNothing = function() list(
  'forceRedoCount'=F,
  'forceRedoNormalCount'=F,
  'forceRedoFit'=F,
  'forceRedoVolcanoes'=F,
  'forceRedoDifferentRegions'=F,
  'forceRedoSNPs'=F,
  'forceRedoVariants'=F,
  'forceRedoNormalSNPs'=F,
  'forceRedoNormalVariants'=F,
  'forceRedoMatchFlag'=F,
  'forceRedoScatters'=F,
  'forceRedoOutputSomatic'=F,
  'forceRedoNewVariants'=F,
  'forceRedoSNPprogression'=F,
  'forceRedoCNV'=F,
  'forceRedoCNVplots'=F,
  'forceRedoSummary'=F,
  'forceRedoStories'=F,
  'forceRedoRiver'=F,
  'forceRedoVEP'=F)

#' returns input that uses saved data if present.
#'
#' @details This function returns the input 'forceRedo' for analyse(), so that all the
#'          analysis is redone, overwriting previously saved data exists.
#' @export
#' @examples
#' \dontrun{
#' cpus=4
#' metaDataFile = '/absolute/path/to/myAnalysis/metaData.txt'
#' vcfFiles = list.files('/absolute/path/to/myAnalysis/vcf', pattern='chr[0-9XYMT]+.vcf$', full.names=T)
#' captureRegionsFile = '/absolute/path/to/myAnalysis/captureRegions.gc.bed'
#' dbSNPdirectory = '/absolute/path/to/myAnalysis/dbSNP'
#' normalDirectory = '/absolute/path/to/myAnalysis/normal'
#' Rdirectory = '/absolute/path/to/myAnalysis/R/'
#' plotDirectory = '/absolute/path/to/myAnalysis/plots/'
#' 
#' #The base quality phred offset. This can be read from fastqc analysis for example.
#' BQoffset = 33
#' genome = 'hg19'
#'
#' inputFiles = list('metaDataFile'=metaDataFile, 'vcfFiles'=vcfFiles, 'normalDirectory'=normalDirectory,
#' 'captureRegionsFile'=captureRegionsFile, 'dbSNPdirectory'=dbSNPdirectory)
#' outputDirectories = list('Rdirectory'=Rdirectory, 'plotDirectory'=plotDirectory)
#' runtimeSettings = list('cpus'=cpus, 'outputToTerminalAsWell'=T)
#'
#' forceRedo = forceRedoEverything()
#' parameters = list('systematicVariance'=0.03, 'maxCov'=150)
#'
#' data = analyse(inputFiles, outputDirectories, settings, forceRedo, runtimeSettings, parameters=parameters)
#' }
forceRedoEverything = function() list(
  'forceRedoCount'=T,
  'forceRedoNormalCount'=T,
  'forceRedoFit'=T,
  'forceRedoVolcanoes'=T,
  'forceRedoDifferentRegions'=T,
  'forceRedoSNPs'=T,
  'forceRedoVariants'=T,
  'forceRedoNormalSNPs'=T,
  'forceRedoNormalVariants'=T,
  'forceRedoMatchFlag'=T,
  'forceRedoScatters'=T,
  'forceRedoOutputSomatic'=T,
  'forceRedoNewVariants'=T,
  'forceRedoSNPprogression'=T,
  'forceRedoCNV'=T,
  'forceRedoCNVplots'=T,
  'forceRedoSummary'=T,
  'forceRedoStories'=T,
  'forceRedoRiver'=T,
  'forceRedoRiver'=T)

#' Wrapper for getting settings, containing defaults for missing values.
#'
#' @param settings A named list of settings.
#' @param name The setting to be retrieved.
#'
#' @details If the setting is present in the list, just return settings[[name]].
#'          Otherwise return a default value.
getSettings = function(settings, name) {
  #if set, return the set value
  if ( name %in% names(settings) ) return(settings[[name]])

  #if not set, return defaults
  if ( name == 'GCcorrection' ) return(TRUE)
  if ( name == 'GCregionCorrection' ) return(FALSE)
  if ( name == 'MAcorrection' ) return(TRUE)
  if ( name == 'sexCorrection' ) return(TRUE)
  if ( name == 'MAFstat' ) return(TRUE)
  if ( name == 'filterFalseSNPregions' ) return(TRUE)
  if ( name == 'useMAFforRegionFinding' ) return(TRUE)
  if ( name == 'varianceBoost' ) return(TRUE)

  #if no default, return NA with a warning
  warning(paste0('setting ', name, ' not set, and has no default.'))
  return(NA)
}


explainSuperMetaData = function() {
  cat('Please provide a (absolute or relative from current directory) path to a metaData file.\n',
      'This is a tab separated file with the mandatory columns\n',
      'BAM\tVCF\tNAME\tINDIVIDUAL\tNORMAL\tTIMEPOINT\n\n',
      'Each row corresponds to one sample to be analysed.\n\n',
      'BAM: path (absolute or relative from metaData file) to bam file of sample.\n',
      'VCF: path (absolute or relative from metaData file) to vcf file of sample. This should include both somatic and germline variants. The variants will undergo further filtering, so it is not a problem if the .vcf includes false positives, although large number of entries in the .vcf increases run time.\n',
      'NAME: unique identifier of the sample. The names will be converted to standard R names, which involves replacing separators (such as space, dash or underscore) with dots.\n',
      'INDIVIDUAL: unique identifier of the indiviudal the sample comes from. This information is used pair up matched normal samples for somatic mutation filtering, and for identifying germline heterozygous SNPs. This column is also used to group samples that should be compared to each other: all samples from the same individual are analysed together, and mutations are tracked between the samples.\n',
      'NORMAL: Should be YES if the sample is normal, NO otherwise. A normal sample is assumed to have no somatic mutations. In practice, a tumor burden below a few percent works effectively as a normal, but tumor burdens above 10% in a sample marked as normal can cause severely reduced sensitivity in cancer samples from the same indiviudal.\n',
      'TIMEPOINT: Mostly used as label in plots together with the sample. Can be Diagnosis, Relapse, resistant or similar labels. Does not influence the analyss otherwise.')  #Probably time to remove the TIMEPOINT column from required...
}

explainSuperCaptureRegions = function() {
  cat('Please provide a (absolute or relative from current directory) path to a captureRegions file.\n',
      'This should be a .bed file: a tab separated file with the first three columns chromosome, start, end.\n')
}

explainSuperNormalDirectory = function() {
  cat('Please provide a (absolute or relative from current directory) path to a normalDirectory.\n',
      'The directory should contain (links to) .bam and .bam.bai for the pool of normal samples.\n',
      'These samples can, but dont have to, be matched normals of the analysed samples.\n',
      'They do need to be from the same capture though, and better if they share error sources with the cancer samples.\n',
      'At least two are required (to estimate variance), but more are better (although slower).\n')
}

explainSuperDbSNPdirectory = function() {
  cat('Please provide a (absolute or relative from current directory) path to the dbSNP directory.\n',
      'This directory is provided with the example for hg19.\n')
}

explainSuperReference = function() {
  cat('Please provide a (absolute or relative from current directory) path to the reference.\n',
      'This should be a fasta file, with an index.\n')
}

#' sets up and checks the input files for superFreq
#'
#' @details Returns the input files in a format to be passed to analyse(). Returns detailed information about the required input files if missing.
#' @export
#' @examples
#' superInputFiles()
superInputFiles = function(metaData='', captureRegions='', normalDirectory='', dbSNPdirectory='', reference='', normalCoverageDirectory='') {
  if ( metaData == '' )
    explainSuperMetaData()
  
  if ( captureRegions == '' )
    explainSuperCaptureRegions()
  
  if ( normalDirectory == '' )
    explainSuperNormalDirectory()
  
  if ( dbSNPdirectory == '' )
    explainSuperDbSNPdirectory()
  
  if ( reference == '' )
    explainSuperReference()
  
  if ( normalCoverageDirectory == '' )
    normalCoverageDirectory = normalDirectory

  ret = list('metaDataFile'=normalizePath(metaData), 'normalDirectory'=normalizePath(normalDirectory),
          'normalCoverageDirectory'=normalizePath(normalCoverageDirectory), 'reference'=normalizePath(reference),
          'captureRegionsFile'=normalizePath(captureRegions), 'dbSNPdirectory'=normalizePath(dbSNPdirectory))

  if ( any(ret=='') ) {
    error('There are required input files not provided.\n')
    return()
  }

  return(ret)
}


#' Sets up and checks the output directories. 
#'
#' @details Outputs some infomration about what the paths should be if called empty.
#' @export
#' @examples
#' superOutputDirectories()
superOutputDirectories = function(Rdirectory='', plotDirectory='') {

  if ( Rdirectory == '' )
    cat('Please provide a (absolute or relative from current directory) path to save results in .Rdata format.\n')

  if ( plotDirectory == '' )
    cat('Please provide a (absolute or relative from current directory) path to save plots in.\n')

  
  parentRdirectory = dirname(Rdirectory)
  if ( !file.exists(parentRdirectory) ) stop('Parent directory of the R directory doesnt exist. Aborting.')
  if ( !file.exists(Rdirectory) ) dir.create(Rdirectory)
  parentPlotDirectory = dirname(plotDirectory)
  if ( !file.exists(parentPlotDirectory) ) stop('Parent directory of the plot directory doesnt exist. Aborting.')
  if ( !file.exists(plotDirectory) ) dir.create(plotDirectory)
  return(list('Rdirectory'=normalizePath(Rdirectory), 'plotDirectory'=normalizePath(plotDirectory)))
}


#' checks if a file exists, and creates an error if it doesnt. 
#'
#' @details Internal method to give approapriate errors for missing files.
#' @examples
#' requireFileExists('a/path/to/a/file.txt')
requireFileExists = function(file) {
  if ( !file.exists(file) ) stop('required file', file, 'not found.\n')
  return(TRUE)
}




#' Returns default settings 
#'
#' @details feed the output of this function to analyse, or just call superFreq.
#' @examples
#' defaultSuperSettings()
defaultSuperSettings = function(genome='', BQoffset='') {
  if ( genome == '' ) {
    cat('Defaulting genome to hg19.\n')
    genome = 'hg19'
  }
  if ( BQoffset == '' ) {
    cat('Defaulting base quality offset to 33.\n')
    BQoffset = 33
  }

  return(list(genome=genome, BQoffset=BQoffset))
}


#' Returns default runtime settings 
#'
#' @details feed the output of this function to analyse, or just call superFreq.
#' @examples
#' defaultSuperRuntimeSettings()
defaultSuperRuntimeSettings = function(cpus='', outputToTerminalAsWell='') {
  if ( cpus == '' ) {
    cat('Defaulting maximum parallel cpus used to 3.\n')
    cpus = 4
  }
  if ( outputToTerminalAsWell == '' ) {
    cat('Defaulting outputToTerminalAsWell to TRUE.\n')
    outputToTerminalAsWell = TRUE
  }

  return(list(cpus=cpus, outputToTerminalAsWell=outputToTerminalAsWell))
}


#' Returns default parameters 
#'
#' @details feed the output of this function to analyse, or just call superFreq.
#' @examples
#' defaultSuperParameters()
defaultSuperParameters = function(systematicVariance='', maxCov='', cloneDistanceCut='') {
  if ( systematicVariance == '' ) {
    cat('Defaulting systematicVariance used to 0.03.\n')
    systematicVariance = 0.03
  }
  if ( maxCov == '' ) {
    cat('Defaulting maxCov to 150.\n')
    maxCov = 150
  }
  if ( cloneDistanceCut == '' ) {
    cat('Defaulting cloneDistanceCut to 2.3.\n')
    cloneDistanceCut = -qnorm(0.01)
  }

  return(list(systematicVariance=systematicVariance, maxCov=maxCov, cloneDistanceCut=cloneDistanceCut))
}


ensureDirectoryExists = function(dir) {
  if ( !file.exists(dir) ) dir.create(dir)
}


#' Downloads resources and template for a standard analysis
#'
#' @param analysisDirectory The directory where the analysis will be made. The directory will be created
#'                          if needed, but the parent directory has to exist.
#'
#' @details This function downloads dbSNP and COSMIC resources for hg19, condensed to what superFreq needs.
#'          It also sets up template files and directories to help set up the project specific data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  
#' downloadTemplate('~/mySuperFreqAnalaysis')
#'
#' }
downloadTemplate = function(analysisDirectory) {
  analysisDirectory = gsub('\\/$', '', analysisDirectory)
  parentDirectory = dirname(analysisDirectory)
  if ( !dir.exists(parentDirectory) ) stop('Parent directory ', parentDirectory, ' does not exist.')

  ensureDirectoryExists(analysisDirectory)

  cat('Downloading:\n')
  cat('COSMIC..')
  cosmicDirectory = paste0(analysisDirectory, '/COSMIC')
  ensureDirectoryExists(cosmicDirectory)
  download.file('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/cosmicCounts.Rdata',
                paste0(cosmicDirectory, '/cosmicCounts.Rdata'), method='curl')
  download.file('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/allCosmicCounts.Rdata',
                paste0(cosmicDirectory, '/allCosmicCounts.Rdata'), method='curl')
  
  cat('dbSNP..')
  dbDirectory = paste0(analysisDirectory, '/dbSNP')
  ensureDirectoryExists(dbDirectory)
  for ( chr in c(as.character(1:22), 'X', 'Y', 'MT') )
    download.file(paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/ds_flat_ch', chr, '.Rdata'),
                  paste0(dbDirectory, '/ds_flat_ch', chr, '.Rdata'), method='curl')


  
}
