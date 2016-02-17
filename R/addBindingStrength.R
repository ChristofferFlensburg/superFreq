

#' Adds binding strength (essentailly same as GC content) to a bedfile.
#'
#' @param bedFile The bedfile.
#' @param fastaFile The fasta file.
#'
#' @details Adds the bidning strength as a column to the bedfile, saving it as a copy, replacing .bed with .dn.bed
addBindingStrength = function(bedFile, fastaFile, genome='hg19', forceRedo=F) {
  outFile = gsub('.bed$', '.dn.bed', bedFile)
  if ( !forceRedo & file.exists(outFile) ) {
    cat('Output file', outFile, 'already exists. Run with forceRedo=TRUE to overwrite.\n')
    return()
  }


  catLog('Loading genome from', fastaFile)
  reference = importReferenceGenome(file=fastaFile)
  names(reference) = gsub('chr', '', names(reference))
  catLog(' done.\n')

  catLog('Loading capture regions from', bedFile)
  cr = read.table(bedFile, fill=T, quote='')
  chrL = chrLengths(genome)
  use = gsub('chr', '',cr$V1) %in% names(chrL)
  cr = cr[use,]
  cr = GRanges(ranges = IRanges(start=cr$V2, end = cr$V3),
    seqnames=gsub('chr', '',cr$V1), seqlengths=chrL, region = gsub(',', '', gsub("\"", '', as.character(cr$V4))))
  names(cr) = cr$region
  catLog('done.\n')

  catLog('Calculating binding strength over capture regions..')
  scores = rangesToBindingScore(as.character(seqnames(cr)), start(cr), end(cr), reference)
  outDF = data.frame(seqnames(cr), start(cr), end(cr), names(cr), scores$baseScores, scores$dibaseScores)
  catLog('done.\n')

  catLog('Printing to', outFile)
  write.table(outDF, file=outFile, quote=F, col.names=F, row.names=F)
  catLog(' done!\n')
}







dibaseStrength = function() {
  return(c('AA'=5, 'AT'=7, 'AC'=10, 'AG'=8,
           'TA'=4, 'TT'=5, 'TC'=8, 'TG'=7,
           'CA'=7, 'CT'=8, 'CC'=11, 'CG'=10,
           'GA'=8, 'GT'=10, 'GC'=13, 'GG'=11))
}
baseStrength = function() {
  return(c('A'=0, 'T'=0, 'C'=1, 'G'=1))
}

dibaseScore = function(seq) {
  if ( nchar(seq) == 0 ) return(-1)
  bM = dibaseStrength()
  diNucs = substring(seq, 1:(nchar(seq)-1), 2:nchar(seq))
  return(mean(bM[diNucs]))
}
baseScore = function(seq) {
  if ( nchar(seq) == 0 ) return(-1)
  bM = baseStrength()
  nucs = substring(seq, 1:nchar(seq), 1:nchar(seq))
  return(mean(bM[nucs]))
}


importReferenceGenome = function(file='~/data/reference/hg19/hg19.fa') {
  ret = read.fasta(file, seqtype='DNA')
  return(ret)
}


rangesToBindingScore = function(chr, start, end, reference) {
  bS = c()
  dbS = c()
  for ( ch in unique(chr) ) {
    cat(ch, '..', sep='')
    chStart = start[chr==ch]
    chEnd = end[chr==ch]
    ranges = lapply(1:length(chStart), function(i) c(chStart[i], chEnd[i]))
    cat('got ', length(ranges), ' ranges.', sep='')
    cat('reference chr..', sep='')
    chRef = reference[[ch]]
    cat('extracting sequences..', sep='')
    chSeqs = sapply(ranges, function(range) toupper(do.call(paste0, as.list(chRef[(range[1]-100):(range[2]+100)]))))
    cat('trimming..', sep='')
    chSeqs = gsub('N', '', chSeqs)
    cat('naming..', sep='')
    names(chSeqs) = paste0(ch, ':', chStart, '-', chEnd)
    cat('base scores..', sep='')
    bS = c(bS, sapply(chSeqs, baseScore))
    cat('dibase scores..', sep='')
    dbS = c(dbS, sapply(chSeqs, dibaseScore))
    cat('done.\n')
  }
  return(list('baseScores'=bS, 'dibaseScores'=dbS))
}
