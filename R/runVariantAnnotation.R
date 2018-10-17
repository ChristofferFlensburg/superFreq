

#' Run VariantAnnotation on the provided samples
#'
#' @param qs list of data.frame: The variants to be annotated by sample
#' @param genome character: hg19, hg38 or mm10.
#' @param cpus integer: The number of cpus to be used as most.
#' @param annotationDirectory character: Where the annotation resources are taken from.
#' @param reference character: Path to the reference fasta.
#'
#' @details This function uses the VariantAnnotation package to annotate somatic variants.
#'          PolyPHen, SIFT and exon number are removed when switching from VEP.
annotateSomaticQs = function(qs, genome='hg19', resourceDirectory='superFreqResources', reference, cpus=1) {

  catLog('Running VariantAnnotation.\n')
  #extract and sort unique variants across samples
  somQlist = lapply(qs, function(q) q[q$somaticP > 0,])
  somQ = do.call(rbind, somQlist)
  if ( any(is.na(somQ$x) | is.na(somQ$somaticP)) ) stop('NA x or somaticP in annotateSomaticQs.')
  somQ = somQ[!duplicated(paste0(somQ$x, somQ$variant)),]
  rownames(somQ) = paste0(somQ$x, somQ$variant)
  #this is the same ordering as in q, so subsetting rows  of somQ will match back to q without sorting.
  somQ = somQ[order(somQ$x, somQ$variant),]

  #do the actual annotation
  somQ = annotateQ(somQ, genome=genome, resourceDirectory=resourceDirectory, reference=reference, cpus=cpus)

  #put annotation back in to the samples
  annotatedQs = lapply(qs, function(q) {
    #this assumes no duplicated rows. Shouldnt be the case, but who knows.
    if ( any(duplicated(paste0(q$x, q$variant))) ) stop('duplicated q rows in annotateSomaticQs.')
    
    #first subsample annotated variants. This will make them in same order as q[somatics,]
    somatics = q$somaticP > 0
    subSomQ = somQ[rownames(somQ) %in% rownames(q)[somatics],]
    if ( !all(rownames(q)[somatics] == rownames(subSomQ)) ) stop('q and subSomQ not matching rownames in annotateSomaticQs')
    
    #set up blank data for the variants not annotated.
    q = addVAnullAnnotation(q, genome)
    
    #fill in annotated variants
    q[somatics, annotationColumns(genome)] = subSomQ[,annotationColumns(genome)]

    return(q)
  })

  return(annotatedQs)
}


#' Run VariantAnnotation on the provided variants
#'
#' @param q data.frame: The variants to be annotated.
#' @param genome character: hg19, hg38 or mm10.
#' @param cpus integer: The number of cpus to be used as most.
#' @param annotationDirectory character: Where the annotation resources are taken from.
#' @param reference character: Path to the reference fasta.
#'
#' @details This function calls VEP on the output from outputSomaticVariants. For this, VEP needs to be callable by system('vep').
#'
#' @importFrom GenomicFeatures makeTxDb
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom VariantAnnotation locateVariants
#' @importFrom VariantAnnotation CodingVariants
#' @importFrom VariantAnnotation AllVariants
#' @importFrom VariantAnnotation predictCoding
#' @importFrom Rsamtools FaFile
#' @importFrom VariantAnnotation locateVariants
#' @importFrom VariantAnnotation VCF
#' @importFrom SummarizedExperiment rowRanges
annotateQ = function(q, genome='hg19', resourceDirectory='superFreqResources', reference, cpus=1) {

  #if empty, add columns and return.
  if ( nrow(q) == 0 ) {
    q = addVAnullAnnotation(q, genome)
    return(q)
  }
  
  #load relevant annotation data
  txdb = ''
  dump = ''
  catLog('Loading annotation dump...')
  if ( genome == 'hg19' ) {
    load(paste0(resourceDirectory, '/annotation/dumpHg19.Rdata'))
    dump = dumpHg19
  }
  if ( genome == 'hg38' ) {
    load(paste0(resourceDirectory, '/annotation/dumpHg38.Rdata'))
    dump = dumpHg38
  }
  if ( genome == 'mm10' ) {
    load(paste0(resourceDirectory, '/annotation/dumpMm10.Rdata'))
    dump = dumpMm10
  }
  catLog('done.\n')
  
  #there is overehead for each batch which favours large batches
  #but also risk for memory crashes if batches are too big.
  #so make one batch per cpus, as long as each batch isnt more than 30k variants.
  #this should be safe in terms of memory with 10Gb per cpu. I think. I hope.
  batches = pmax(1, round(nrow(q)/10e3))
  breaks = round((0:batches)/batches*(nrow(q)) + 1)
  qList = lapply(1:batches, function(i) q[breaks[i]:(breaks[i+1]-1),])
  catLog('Splitting up ',  nrow(q), ' variants for parallelisation into ', batches, ' batches.\n', sep='')

  #set up link to fasta for this thread, and check if "chr" is present in chr names.
  catLog('Setting up data bases..')
  fafile = FaFile(reference)
  hasChr = grepl('^chr', as.character(seqnames(scanFaIndex(fafile))[1]))
  if ( !hasChr ) dump = removeChrFromDump(dump)
  
  #create data base for this thread
  txdb = GenomicFeatures::makeTxDb(transcripts=dump$transcripts, splicings=dump$splicings, genes=dump$genes, chrominfo=dump$chrominfo)
  catLog('done.\n')
  
  catLog('Running annotation by batch')
  qList = lapply(qList, function(q) {
    catLog('.')
        
    #get the variants into the vcf format. Not sure it handles indels correctly...
    vcf = superFreq:::qToVCF(q, genome=genome, seqlengths=seqlengths(txdb), addChrToSeqnames=hasChr, fafile=fafile)
    rd = SummarizedExperiment::rowRanges(vcf)
    
    #these throw a warning about out of bound granges, not sure why, seems to return ok values.
    suppressWarnings(suppressMessages({allvar = VariantAnnotation::locateVariants(rd, txdb, VariantAnnotation::AllVariants())}))
    suppressWarnings({coding = VariantAnnotation::predictCoding(vcf, txdb, seqSource=fafile)})
    
    #rank and select most severe hit for each variant
    coding$SEVERITY = superFreq:::VAconsequenceToSeverityRank(as.character(coding$CONSEQUENCE))
    allvar$SEVERITY = superFreq:::VAconsequenceToSeverityRank(as.character(allvar$LOCATION))
    
    #select most severe and merge to q
    q = superFreq:::addMostSevereHit(q, allvar, coding, genome)
    return(q)
  })
  catLog('done.\n')

  CCGDsummary = ''
  if ( genome %in% c('mm10') ) {
    catLog('Preprocessing CCGD data...')
    CCGDsummary = preprocessCCGD(resourceDirectory=resourceDirectory)
    catLog('done.\n')
  }
  
  catLog('Matching to known cancer variants by batch')
  qList = mclapply(qList, function(q) {
    catLog('.')
    #add frequent cancer genes from COSMIC or CCGD depending on human or mouse
    if ( genome %in% c('hg19', 'hg38') )
      q = addCOSMICannotation(q, genome=genome, resourceDirectory=resourceDirectory)
    if ( genome %in% c('mm10') )
      q = superFreq:::addCCGDannotation(q, CCGDsummary)

    #add clinvar annotation
    if ( genome %in% c('hg19', 'hg38') )
      q = addClinvarAnnotation(q, genome=genome, resourceDirectory=resourceDirectory)

    return(q)
  }, mc.cores=cpus)
  q = do.call(rbind, qList)
  catLog('done.\n')

  return(q)
}

addMostSevereHit = function(q, allvar, coding, genome) {
  #add/clear columns
  q = addVAnullAnnotation(q, genome)

  #sort by severity and remove duplicates, then sort on q ID again.
  allvar = allvar[order(as.numeric(allvar$SEVERITY))]
  allvar = allvar[!duplicated(allvar$QUERYID)]
  allvar = allvar[order(allvar$QUERYID)]
  coding = coding[order(as.numeric(coding$SEVERITY))]
  coding = coding[!duplicated(coding$QUERYID)]
  coding = coding[order(coding$QUERYID)]

  #set up data frame with relevant columns
  allvarDF =
    data.frame(severity=as.numeric(allvar$SEVERITY),
               type=as.character(allvar$LOCATION), stringsAsFactors=F)
  codingDF =
    data.frame(severity=as.numeric(coding$SEVERITY),
               type=as.character(coding$CONSEQUENCE),
               AApos=sapply(coding$PROTEINLOC, function(x) as.numeric(x)[1]),
               AAbefore=sapply(coding$REFAA, function(x) as.character(x)),
               AAafter=sapply(coding$VARAA, function(x) as.character(x)), stringsAsFactors=F)
  
  
  #insert into q
  q[allvar$QUERYID, colnames(allvarDF)] = allvarDF
  q[coding$QUERYID, colnames(codingDF)] = codingDF

  return(q)
}

#' Transforms a variant q data frame into granges
#'
#' @importFrom GenomicRanges GRanges
qToGRanges = function(q, genome='hg19', addChrToSeqnames=F, seqlengths) {
    chr = xToChr(q$x, genome)
    if ( addChrToSeqnames ) chr = paste0('chr', chr)
    pos = xToPos(q$x, genome)
    gr = GRanges(seqnames=chr, ranges = IRanges(start=pos, end=pos+nchar(q$reference)-1), strand='*', seqlengths=seqlengths)
    return(gr)
}

qToVCF = function(q, genome, addChrToSeqnames, seqlengths, fafile) {  
  REF = DNAStringSet(q$reference)
  ALT = CharacterList(as.list(q$variant))
  QUAL = rep(40, nrow(q)) #placeholder, only want annotation anyway.
  FILTER = rep('PASS', nrow(q)) #placeholder, only want annotation anyway.
  refvarstartendList = convertVariantsToVCF(list(REF, ALT, q$x, q$x))
  REF = refvarstartendList[[1]]
  ALT = refvarstartendList[[2]]

  colData = DataFrame(Samples='sample')

  q$reference = as.character(REF)
  gr = qToGRanges(q, genome=genome, addChrToSeqnames=addChrToSeqnames, seqlengths=seqlengths)
  REF = getSeq(fafile, gr)
  ALT[nchar(q$reference)>1] = substr(as.character(REF[nchar(q$reference)>1]), 1, 1)

  fixed = DataFrame(REF=REF, ALT=ALT, QUAL=QUAL, FILTER=FILTER)
  
  vcf = VariantAnnotation::VCF(rowRanges = gr, fixed=fixed, colData=colData, collapsed=TRUE, verbose=F)

  return(vcf)
}



#this is in the outputsomaticVariants.R, so can be removed when merged into the package
convertVariantsToVCF = function(refvarstartendList) {
  reference = refvarstartendList[[1]]
  variant = refvarstartendList[[2]]
  start = refvarstartendList[[3]]
  end = refvarstartendList[[4]]
  if ( any(grepl('-', variant)) ) {
    deletions = grepl('-', variant)
    nDel = as.numeric(gsub('-', '', variant[deletions]))
    reference[deletions] = sapply(nDel, function(n) do.call(paste0, as.list(rep('N', n+1))))
    variant[deletions] = 'N'
    start[deletions] = start[deletions]+1
    end[deletions] = end[deletions]+nDel
  }
  if ( any(grepl('\\+', variant)) ) {
    insertions = grepl('\\+', variant)
    nIns = nchar(as.character(variant[insertions]))-1
    variant[insertions] = paste0(substr(reference[insertions], 1, 1), gsub('\\+', '', variant[insertions]))
  }

  return(list(reference, variant, start, end))
}

#converts the consequence to a severity rank.
VAconsequenceToSeverityRank = function(consequence) {
  #'coding' is used by the location run to cover all the coding variants, and will be replaced
  #by the consequence from the coding run, which is why it's ranked last. Should never be top hit.
  #these number are set to match the VEP ranking.
  #some of these categories include multiple VEP categories. It's not 1 to 1, so some shifting is done.
  #11 still separates coding from non-coding variants.
  rank = c(
    'spliceSite'=2,
    'nonsense'=4,
    'frameshift'=5,
    'nonsynonymous'=10,
    'synonymous'=16,
    'promoter'=19,
    'threeUTR'=20,
    'fiveUTR'=21,
    'intron'=22,
    'intergenic'=30,
    'coding'=100)
  ret = rank[consequence]
  ret[is.na(ret)] = 100
  return(ret)
}

#converts the consequence to a severity rank.
VAseverityToConsequence = function(severity) {
  #'coding' is used by the location run to cover all the coding variants, and will be replaced
  #by the consequence from the coding run, which is why it's ranked last. Should never be top hit.
  #these number are set to match the VEP ranking.
  #some of these categories include multiple VEP categories. It's not 1 to 1, so some shifting is done.
  #11 still separates coding from non-coding variants.
  consequence = c(
    '2'='spliceSite',
    '4'='nonsense',
    '5'='frameshift',
    '10'='nonsynonymous',
    '16'='synonymous',
    '19'='promoter',
    '20'='threeUTR',
    '21'='fiveUTR',
    '22'='intron',
    '30'='intergenic',
    '100'='unknown')
  ret = consequence[as.character(severity)]
  ret[is.na(ret)] = 'unknown'
  return(ret)
}



#the expected EXTRA columns added by annotation.
annotationColumns = function(genome) {
  if ( genome %in% c('hg19', 'hg38') )
    return(c('severity', 'type', 'isCosmicCensus', 'AApos', 'AAbefore', 'AAafter', 'cosmicVariantRank', 'cosmicVariantFrequency', 'cosmicGeneRank', 'cosmicGeneFrequency', 'ClinVar_ClinicalSignificance', 'ClinVar_OriginSimple', 'ClinVar_ReviewStatus', 'ClinVar_VariationID'))
  if ( genome %in% c('mm10') )
    return(c('severity', 'type', 'AApos', 'AAbefore', 'AAafter', 'isCCGD', 'CCGDstudies', 'CCGDcancerTypes', 'CCGDcosmic', 'CCGDcgc', 'CCGDranks'))
}


addVAnullAnnotation = function(q, genome) {
  q$severity = rep(100, nrow(q))
  q$type = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$isCosmicCensus = rep(F, nrow(q))
  q$AApos = rep('', nrow(q))
  q$AAbefore = rep('', nrow(q))
  q$AAafter = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$cosmicVariantRank = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$cosmicVariantFrequency = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$cosmicGeneRank = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$cosmicGeneFrequency = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$ClinVar_ClinicalSignificance = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$ClinVar_OriginSimple = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$ClinVar_ReviewStatus = rep('', nrow(q))
  if ( genome %in% c('hg19', 'hg38') ) q$ClinVar_VariationID = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$isCCGD = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$CCGDstudies = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$CCGDcancerTypes = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$CCGDcosmic = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$CCGDcgc = rep('', nrow(q))
  if ( genome %in% c('mm10') ) q$CCGDranks = rep('', nrow(q))
  return(q)
}

#this function adds COSMIC annotation to q
addCOSMICannotation = function(q, genome, resourceDirectory) {
  #first load and save the census gene names
  countsFile = paste0(resourceDirectory, "/COSMIC/cosmicCounts.Rdata")
  load(countsFile)
  censusGenes = names(cosmicCounts$geneCounts)
  q$isCosmicCensus = q$inGene %in% censusGenes

  #then load in the full COSMIC information
  countsFile = paste0(resourceDirectory, "/COSMIC/allCosmicCounts_", genome, ".Rdata")
  load(countsFile)

  #get the rank by gene, and filter the transcript-specific [SYMBOL]_[ENSEMBL] entries.
  genesDF =
    data.frame(gene=names(cosmicCounts$geneDensity),
               density=as.numeric(cosmicCounts$geneDensity), stringsAsFactors=F)
  genesDF = genesDF[!grepl('_ENST', genesDF$gene),]
  genesDF$rank = rank(-genesDF$density)

  #link gene rates to q
  genesDF = genesDF[genesDF$gene %in% q$inGene,] #this speeds up the random access next lines
  rownames(genesDF) = genesDF$gene
  genesDF = genesDF[q$inGene,]
  q$cosmicGeneMPMPB = genesDF$density
  q$cosmicGeneRank = genesDF$rank

  #get the rate and rank by position
  xr = cosmicCounts$xRates
  xr$rank = rank(-xr$rates)

  #link position rates to q
  xr = xr[xr$x %in% q$x,] #this speeds up the random access next lines
  rownames(xr) = xr$x
  xr = xr[as.character(q$x),]
  q$cosmicVariantFrequency = xr$rates
  q$cosmicVariantRank = xr$rank

  return(q)
}

preprocessCCGD = function(resourceDirectory) {
  CCGDdata = read.table(paste0(resourceDirectory, '/COSMIC/CCGD_export.csv'), sep=',', header=T, stringsAsFactors=F)
  censusGenes = unique(CCGDdata$Mouse.Symbol)
  
  CCGDsummary = sapply(censusGenes, function(gene) {
    subData = CCGDdata[CCGDdata$Mouse.Symbol == gene,]
    studies = subData$Studies[1]
    humanGene = subData$Human.Symbol[1]
    
    cancerTypes = sort(table(subData$Cancer.Type), decreasing=T)
    names(cancerTypes) = gsub(' Cancer', '',names(cancerTypes))
    cancerTypes = gsub(';$','',do.call(paste0, as.list(paste0(names(cancerTypes), ":", cancerTypes, ';'))))
    
    isInCosmic = subData$COSMIC[1] == 'Yes'
    isInCGC = subData$CGC[1] == 'Yes'
    
    ranks = table(c('A', 'B', 'C', 'D', subData$Relative.Rank))
    ranks[c('A', 'B', 'C', 'D')] = ranks[c('A', 'B', 'C', 'D')] - 1
    names(ranks) = gsub(' Cancer', '',names(ranks))
    ranks = gsub(';$','',do.call(paste0, as.list(paste0(names(ranks), ":", ranks, ';'))))
    
    return(c('studies'=studies, 'cancerTypes'=cancerTypes, 'cosmic'=isInCosmic,
             'cgc'=isInCGC, 'ranks'=ranks, 'humanGene'=humanGene))
  })
  CCGDsummary = as.data.frame(t(CCGDsummary), stringsAsFactors=F)
  rownames(CCGDsummary) = censusGenes
  CCGDsummary$studies = as.numeric(CCGDsummary$studies)
  CCGDsummary$cosmic = as.logical(CCGDsummary$cosmic)
  CCGDsummary$cgc = as.logical(CCGDsummary$cgc)

  return(CCGDsummary)
}

#adds CCGD annotation to q from the q$inGene column.
addCCGDannotation = function(q, CCGDsummary) {
  
  censusGenes = unique(CCGDdata$Mouse.Symbol)
  isCCGD = q$inGene %in% censusGenes
    
  CCGDstudies = rep(0, nrow(q))
  CCGDcancerTypes = rep('', nrow(q))
  CCGDcosmic = rep('', nrow(q))
  CCGDcgc = rep('', nrow(q))
  CCGDranks = rep('', nrow(q))
  
  CCGDstudies[isCCGD] = CCGDsummary[q$inGene[isCCGD],]$studies
  CCGDcancerTypes[isCCGD] = CCGDsummary[q$inGene[isCCGD],]$cancerTypes
  CCGDcosmic[isCCGD] = ifelse(CCGDsummary[q$inGene[isCCGD],]$cosmic, 'YES!', 'no')
  CCGDcgc[isCCGD] = ifelse(CCGDsummary[q$inGene[isCCGD],]$cgc, 'YES!', 'no')
  CCGDranks[isCCGD] = CCGDsummary[q$inGene[isCCGD],]$ranks

  q$isCCGD = isCCGD
  q$CCGDstudies = CCGDstudies
  q$CCGDcancerTypes = CCGDcancerTypes
  q$CCGDcosmic = CCGDcosmic
  q$CCGDcgc = CCGDcgc
  q$CCGDranks = CCGDranks
  
  return(q)
}


addClinvarAnnotation = function(q, genome, resourceDirectory) {
  if ( genome == 'hg19' ) load(paste0(resourceDirectory, '/COSMIC/clinvar19.Rdata'))
  if ( genome == 'hg38' ) load(paste0(resourceDirectory, '/COSMIC/clinvar38.Rdata'))

  clinvar = clinvar[rownames(clinvar) %in% rownames(q),]
  colnames(clinvar) = paste0('ClinVar_', colnames(clinvar))
  q[,colnames(clinvar)] = ''
  q[rownames(clinvar),colnames(clinvar)] = clinvar

  return(q)
}


removeChrFromDump = function(dump) {
  dump$transcripts$tx_chrom = sub('chr', '', dump$transcripts$tx_chrom)
  dump$splicings$exon_chrom = sub('chr', '', dump$splicings$exon_chrom)
  dump$chrominfo$chrom = sub('chr', '', dump$chrominfo$chrom)
  return(dump)
}


#this is the code used to set up the .Rdata objects that are stored on WEHI. Not called by superFreq.
#all the resources (dbSNP, ExAC, COSMIC, CCGD, standard capture regions, annotation)
#are preprocessed and stored at WEHI for speed and to rely on a single online connection only.
#so this code is only a FYI for how the WEHI data is produced,
#and for maintainers to update or regenerate the WEHI resources.
createOfflineTxDbs = function(annotationDirectory) {
  
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  x = TxDb.Hsapiens.UCSC.hg19.knownGene
  dumpHg19 = as.list(x)
  save(dumpHg19, file=paste0(annotationDirectory, '/dumpHg19.Rdata'))

  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  x = TxDb.Hsapiens.UCSC.hg38.knownGene
  dumpHg38 = as.list(x)
  save(dumpHg38, file=paste0(annotationDirectory, '/dumpHg38.Rdata'))


  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  x = TxDb.Mmusculus.UCSC.mm10.knownGene
  dumpMm10 = as.list(x)
  save(dumpMm10, file=paste0(annotationDirectory, '/dumpMm10.Rdata'))

}

#this function preprocesses the tab separated variant summary table from clinvar.
#this is not done by superFreq, but preprocessed and downloaded with other resources.
preprocessClinvar = function(clinvarFile, saveDir) {
  data = read.table('~/majewski/resources/clinvar/variant_summary.txt', sep='\t', header=T,
    fill=F, stringsAsFactors=F, comment.char='', quote='')

  #split by hg19 and hg38
  dat19 = data[data$Assembly == 'GRCh37',]
  dat38 = data[data$Assembly == 'GRCh38',]

  #get x and variant and rowname
  dats = lapply(list('hg19'=dat19, 'hg38'=dat38), function(dat) {
    genome = 'hg19'
    if ( dat$Assembly[1] == 'GRCh38' ) genome = 'hg38'

    dat = dat[dat$Chromosome %in% names(chrLengths(genome)),]
    dat = dat[dat$Type %in% c('deletion', 'insertion', 'single nucleotide variant'),]
    dat = dat[dat$AlternateAllele != 'na',]
    
    x = chrToX(dat$Chromosome, dat$Start, genome)
    variant = dat$AlternateAllele
    variant = ifelse(dat$Type == 'deletion', paste0('-', nchar(dat$ReferenceAllele)), variant)
    variant = ifelse(dat$Type %in% c('insertion'), paste0('+', dat$AlternateAllele), variant)
    varname = paste0(x, variant)

    ret = dat[,c('ClinicalSignificance', 'OriginSimple', 'ReviewStatus', 'VariationID')]
    ret = ret[!duplicated(varname),]
    rownames(ret) = varname[!duplicated(varname)]

    return(ret)
    })
  
  clinvar = dats$hg19
  save(clinvar, file=paste0(saveDir, '/clinvar19.Rdata'))
  clinvar = dats$hg38
  save(clinvar, file=paste0(saveDir, '/clinvar38.Rdata'))
}


#set up for example during development
if ( F ) {
  #uval melanoma from TCGA
  setwd('~/gdc/TCGA/analysis/TCGA-YZ-A985/')
  library(superFreq)
  source('~/home/Rscripts/CNVcaller/runVariantAnnotation.R')
  library(VariantAnnotation)
  genome='hg38'
  resourceDirectory = '~/temp/sfRes'
  reference = '~/home/leukemia_genomics/data/reference/hg38/hg38.fa'
  cpus = 5

  load('~/gdc/TCGA/analysis/TCGA-YZ-A985/WXS_R/allVariants.Rdata')
  q = allVariants$variants$variants$TCGA.YZ.A985.PrimaryTumor.WXS.194
  q = q[q$somaticP > 0,]
  qs = allVariants$variants$variants
  catLog = cat


  #an internal dlbcl sample
  library(superFreq)
  source('~/home/Rscripts/CNVcaller/runVariantAnnotation.R')
  library(VariantAnnotation)
  genome='hg19'
  resourceDirectory = '~/gdc/superFreq/testCNVdlbcl/superFreqResources'
  reference = '~/gdc/superFreq/testCNVdlbcl/reference/hg19.fa'
  cpus = 5
  load('~/gdc/superFreq/testCNVdlbcl/R/allVariantsPreVEP.Rdata')
  qs = allVariants$variants$variants
  catLog = cat


  #an internal dlbcl sample
  library(superFreq)
  source('~/home/Rscripts/CNVcaller/runVariantAnnotation.R')
  library(VariantAnnotation)
  genome='hg19'
  resourceDirectory = '~/gdc/superFreq/testCNV108/superFreqResources'
  reference = '~/gdc/superFreq/testCNV108/reference/hg19.fa'
  cpus = 5
  load('~/gdc/superFreq/testCNV108/R/patient1/allVariantsPreVEP.Rdata')
  qs = allVariants$variants$variants
  catLog = cat

  #an internal mm10 sample
  library(superFreq)
  source('~/home/Rscripts/CNVcaller/runVariantAnnotation.R')
  library(VariantAnnotation)
  genome='mm10'
  resourceDirectory = '~/gdc/superFreq/testCNVmm/superFreqResources'
  reference = '~/gdc/superFreq/testCNVmm/reference/mm10.fa'
  cpus = 5

  #try running without this to test performance on 60k variants.
  load('~/gdc/superFreq/testCNVmm/R/mouse/allVariantsPreVEP.Rdata')
  qs = allVariants$variants$variants


  #annotate
  qs = annotateSomaticQs(qs=qs, genome=genome, resourceDirectory=resourceDirectory, reference=reference, cpus=cpus)

  #plot scatter
  source('~/home/Rscripts/CNVcaller/makeScatterPlots.R')
  vafScatter(qs[[2]], qs[[1]], plotFlagged=F, cpus=cpus, genome=genome)

  #river
  load('~/gdc/superFreq/testCNVuveal/R/stories.Rdata')
  stories = stories$stories
  load('~/gdc/superFreq/testCNVuveal/R/allVariants.Rdata')
  patVar = allVariants$variants
  patVar$variants = patVar$variants[colnames(stories[[1]]$consistentClusters$cloneStories$stories)]
  plotRiver(cloneTree=stories[[1]]$consistentClusters$cloneTree, cloneStories=stories[[1]]$consistentClusters$cloneStories, storyList=stories[[1]]$consistentClusters$storyList, allStories=stories[[1]]$allConsistent, variants=patVar, genome='hg38')

  #somaticVariants
  load('~/gdc/superFreq/testCNVuveal/R/allVariants.Rdata')
  variants = allVariants$variants
  somatics = outputSomaticVariants(variants, genome='hg38', '~/gdc/superFreq/testCNVuveal/plots', cpus=cpus, forceRedo=T)

  #function that compares calls between VEP and VariantAnnotator
  compareVEP = function(allVariantsFile, genome='hg19') {
    load(allVariantsFile)
    qs = allVariants$variants$variants
    newQs = annotateSomaticQs(qs=qs, genome=genome, resourceDirectory=resourceDirectory, reference=reference, cpus=cpus)
    use = qs[[1]]$severity < 100
    plotColourScatter(qs[[1]]$severity[use], newQs[[1]]$severity[use])
    segments(0,11,100,11)
    segments(11,0,11,100)
    print(table(paste0(qs[[1]]$severity[use], ' - ', newQs[[1]]$severity[use])))
  }
}





