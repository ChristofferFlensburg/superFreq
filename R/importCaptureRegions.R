
#' Imports capture regions from a bed file with GC information.
#'
#' @param file Bedfile of the capture regions with extra columns for gc and binding strength information.
#' @param genome Which genome the bed file is on. 'hg19' and 'mm10' only supported atm.
#' @param gcColumn the column in the bed file that contains the gc content. Default 5 matches output from the provided script at X.
#' @param dnColumn the column in the bed file that contains the binding strength. Default 6 matches output from the provided script at X.
#'
#' @details Reads in the bed file into GRanges format, with extra columns for GC and binding strength.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom biomaRt useMart getBM
#'
#' @examples
#' captureRegions = importCaptureRegions('path/to/my/captureRegions.bed')
#' plotColourScatter(captureRegions$gc, captureRegions$dn, xlab='GC content', ylab='binding strength')
#'
importCaptureRegions = function(bedFile, reference, Rdirectory, genome='hg19', forceRedoCR=F, forceRedoRdirectorySave=F) {
  saveFile = paste0(Rdirectory, '/captureRegions.Rdata')
  if ( file.exists(saveFile) & !forceRedoCR & !forceRedoRdirectorySave ) {
    catLog('Loading capture regions..')
    load(saveFile)
    catLog('done.\n')
    return(captureRegions)
  }
  if ( !grepl('\\.bed$', bedFile) ) stop('Need .bed capture regions.')
  bsFile = gsub('\\.bed$', '.named.dn.bed', bedFile)
  if ( !file.exists(bsFile) | forceRedoCR ) {
    catLog('Couldnt find a binding strength file at ', bsFile, ' Making one now.\n')
    namedFile = gsub('\\.bed$', '.named.bed', bedFile)
    if ( !file.exists(namedFile) | forceRedoCR ) {
      catLog('Couldnt find a named capture region file at ', namedFile, '. Making one now.\n', sep='')
      nameCaptureRegions(bedFile, namedFile, Rdirectory, genome)
      catLog('Done and saved named file, now proceeding to binding strength.\n')
    }
    if ( reference == '' ) stop('Need a reference entry in the input files to calcuate binding strengths.')
    addBindingStrength(namedFile, reference, genome=genome)
  }
  catLog('Loading capture regions...')
  cR = read.table(bsFile, fill=T, quote='')
  chrL = chrLengths(genome)
  use = gsub('chr', '',cR$V1) %in% names(chrL)
  cR = cR[use,]
  captureRegions = GRanges(ranges = IRanges(start=cR$V2, end = cR$V3),
    seqnames=gsub('chr', '',cR$V1), seqlengths=chrL, region = gsub(',', '', gsub("\"", '', as.character(cR$V4))),
    gc = cR[[5]], dn = cR[[6]])

  names(captureRegions) = captureRegions$region
  save(captureRegions, file=saveFile)
  catLog('done.\n')
  return(captureRegions)
}

#helper functions that returns the lengths of the chromosomes of a genome.
humanAllChrLengths = function() {
lengths = c(249250621, 106433, 547496, 243199373, 198022430, 191154276, 590426, 189789, 191469, 180915260, 171115067, 4622290, 4795371,
4610396, 4683263, 4833398, 4611984, 4928567, 159138663, 182896, 146364022, 38914, 37175, 141213431, 90085, 169874, 187035, 36148,
135534747, 135006516, 40103, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 1680828, 37498, 81310, 174588, 41001,
78077248, 4262, 59128983, 92689, 159169, 63025520, 48129895, 27682, 51304566, 155270560, 59373566, 166566, 186858, 164239, 137718,
172545, 172294, 172149, 161147, 179198, 161802, 155397, 186861, 180455, 179693, 211173, 15008, 128374, 129120, 19913, 43691, 27386,
40652, 45941, 40531, 34474, 41934, 45867, 39939, 33824, 41933, 42152, 43523, 43341, 39929, 36651, 38154, 36422, 39786, 38502, 16571)

names(lengths) = c('1','1_gl000191_random','1_gl000192_random','2','3','4','4_ctg9_hap1','4_gl000193_random','4_gl000194_random',
'5','6','6_apd_hap1','6_cox_hap2','6_dbb_hap3','6_mann_hap4','6_mcf_hap5','6_qbl_hap6','6_ssto_hap7','7','7_gl000195_random',
'8','8_gl000196_random','8_gl000197_random','9','9_gl000198_random','9_gl000199_random','9_gl000200_random','9_gl000201_random',
'10','11','11_gl000202_random','12','13','14','15','16','17','17_ctg5_hap1','17_gl000203_random','17_gl000204_random',
'17_gl000205_random','17_gl000206_random','18','18_gl000207_random','19','19_gl000208_random',
'19_gl000209_random', '20', '21', '21_gl000210_random', '22', 'X', 'Y', 'Un_gl000211', 'Un_gl000212', 'Un_gl000213', 'Un_gl000214',
'Un_gl000215', 'Un_gl000216', 'Un_gl000217', 'Un_gl000218', 'Un_gl000219', 'Un_gl000220', 'Un_gl000221', 'Un_gl000222', 'Un_gl000223',
'Un_gl000224', 'Un_gl000225', 'Un_gl000226', 'Un_gl000227', 'Un_gl000228', 'Un_gl000229', 'Un_gl000230', 'Un_gl000231', 'Un_gl000232',
'Un_gl000233', 'Un_gl000234', 'Un_gl000235', 'Un_gl000236', 'Un_gl000237', 'Un_gl000238', 'Un_gl000239', 'Un_gl000240', 'Un_gl000241',
'Un_gl000242', 'Un_gl000243',  'Un_gl000244', 'Un_gl000245', 'Un_gl000246', 'Un_gl000247', 'Un_gl000248', 'Un_gl000249', 'M')

  return(lengths)
}
humanChrLengths = function() {
  humanAllChrLengths()[as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y', 'M'))]
}
mouseChrLengths = function() {
  lengths = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993,
    122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299)

  names(lengths) = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10' ,'11', '12' ,'13', '14' ,'15', '16' ,'17', '18', '19', 'X', 'Y', 'M')

  return(lengths)
}
chrLengths = function(genome='hg19') {
  if ( genome == 'hg19' ) return(humanChrLengths())
  else if ( genome == 'mm10' ) return(mouseChrLengths())
  else stop('chrLengths doesnt know about genome ', genome, '\n')
}


makeCaptureFromFeatureCountRegions = function(outputBed, genome='hg19') {
  ann = getInBuiltAnnotation('hg19')
  entrezGenes = unique(ann$GeneID)
  genes1 = queryMany(entrezGenes[1:10000], species='human')
  genes2 = queryMany(entrezGenes[10001:20000], species='human')
  genes3 = queryMany(entrezGenes[20001:length(entrezGenes)], species='human')
  genes = rbind(as.data.frame(genes1, stringsAsFactors=F), as.data.frame(genes2, stringsAsFactors=F), as.data.frame(genes3, stringsAsFactors=F))

  symbol = genes$symbol
  symbol[is.na(symbol)] = genes$query[is.na(symbol)]
  names(symbol) = genes$query

  chr = annotationToChr(ann, genome)
  start = ann$Start
  end = ann$End
  x1x2 = annotationToX1X2(ann, genome)
  x1 = xToPos(x1x2$x1, genome)
  x2 = xToPos(x1x2$x2, genome)
  ret = rbind(chr, start, pmax(start, end), symbol[as.character(ann$GeneID)])
  ret = ret[,ret[1,] %in% names(chrLengths(genome))]

  write(ret, file=outputBed, sep='\t', ncolumns=4)
}


padCaptureRegions = function(inFile, outFile=gsub('.bed$', '.padded.bed', inFile), padding=100) {
  if ( file.exists(outFile) ) return()
  cr = importCaptureRegions(inFile, gcColumn=3, dnColumn=3)
  outDF = data.frame(seqnames(cr), pmax(1, start(cr)-padding), pmin(end(cr)+padding), names(cr))
  temp = options()$scipen
  options(scipen = 100)
  write.table(outDF, file=outFile, quote=F, col.names=F, row.names=F)
  options(scipen = temp)
}


nameCaptureRegions = function(inFile, outFile, Rdirectory, genome='hg19') {
  cr = read.table(inFile, fill=T, quote='')
  chr = gsub('chr', '',cr$V1)
  start = cr$V2
  end = cr$V3

  use = chr %in% names(chrLengths(genome))
  if ( sum(!use) > 10 ) warning(paste0(sum(!use), ' rows in the capture region file does not match a chromosome in the genome ', genome, ' will use the remaining ', sum(use), ' rows.'))
  catLog('Removing mitohondrial regions from the analysis.\n')
  use = use & !(chr %in% c('M', 'MT'))
  if ( sum(use) == 0 ) stop('No rows in the capture regions start with a chromosome name. Aborting, sorry.')
  chr = as.character(chr[use])
  start = as.numeric(start[use])
  end = as.numeric(end[use])
  x = chrToX(chr, (start+end)/2, genome=genome)
  gene = xToGeneFromDB(x, genome=genome, saveDirectory=Rdirectory)

  outDF = data.frame(chr, start, end, gene)
  catLog('done.\n')

  catLog('Printing to', outFile)
  write.table(outDF, file=outFile, quote=F, col.names=F, row.names=F)
  catLog(' done!\n')
}


#' find the gene closest to the genomic coordinate x
#'
#' @importFrom GenomicRanges findOverlaps as.data.frame as.factor
xToGene = function(x, captureRegions, genome='hg19') {
  if ( length(x) == 0 ) return(c())

  chr = xToChr(x, genome)
  pos = xToPos(x, genome)
  xGR = GRanges(seqnames=chr, ranges = IRanges(start=pos, end=pos), strand='*', seqlengths=chrLengths(genome))

  gene = rep('', length(x))
  fO = findOverlaps(xGR, captureRegions, maxgap=0)
  OLs = cbind(queryHits(fO), subjectHits(fO))
  OLs = OLs[!duplicated(OLs[,1]),,drop=F]
  gene[OLs[,1]] = captureRegions$region[OLs[,2]]
  fO = findOverlaps(xGR[gene == ''], captureRegions, maxgap=100)
  OLs = cbind(queryHits(fO), subjectHits(fO))
  OLs = OLs[!duplicated(OLs[,1]),,drop=F]
  gene[which(gene == '')[OLs[,1]]] = captureRegions$region[OLs[,2]]
  fO = findOverlaps(xGR[gene == ''], captureRegions, maxgap=200)
  OLs = cbind(queryHits(fO), subjectHits(fO))
  OLs = OLs[!duplicated(OLs[,1]),,drop=F]
  gene[which(gene == '')[OLs[,1]]] = captureRegions$region[OLs[,2]]
  fO = findOverlaps(xGR[gene == ''], captureRegions, maxgap=300)
  OLs = cbind(queryHits(fO), subjectHits(fO))
  OLs = OLs[!duplicated(OLs[,1]),,drop=F]
  gene[which(gene == '')[OLs[,1]]] = captureRegions$region[OLs[,2]]

  if ( any(gene=='') ) warning('Failed to find gene names at x=', x[gene==''][1],' in xToGene.')
  gene[gene==''] = '???'

  return(gene)
}


xToGeneFromDB = function(x, genome='hg19', saveDirectory='', verbose=T) {
  if ( length(x) == 0 ) return(c())

  bm = importEnsemblData(x, saveDirectory, genome, verbose=verbose)
  
  chr = xToChr(x, genome)
  pos = xToPos(x, genome)

  #match variant positions to genes ranges
  if ( verbose ) catLog('matching ranges to positions...')
  if ( genome == 'hg19' ) symbolName = 'hgnc_symbol'
  else if ( genome == 'mm10' ) symbolName = 'mgi_symbol'
  else stop('genome ', genome, ' not supported')
  geneRanges =
    GRanges(seqnames=bm$chromosome_name,
            ranges = IRanges(start=bm$start_position, end=bm$end_position, names=bm[[symbolName]]),
            strand='*', seqlengths=chrLengths(genome))
  variantRanges = GRanges(seqnames=chr, ranges = IRanges(start=pos, end=pos, names=x),
    strand='*', seqlengths=chrLengths(genome))
  #allow hits up to 500bp away.
  fO = findOverlaps(variantRanges, geneRanges, maxgap=500)
  OLs = cbind(queryHits(fO), subjectHits(fO))
  #remove duplicates, leaving only lowest index hit for each variant
  OLs = OLs[order(OLs[,2]),,drop=F]
  OLs = OLs[!duplicated(OLs[,1]),,drop=F]
  
  rownames(OLs) = OLs[,1]

  xIndex = as.character(1:length(x))
  hasHit = xIndex %in% rownames(OLs)
  genes = rep('?', length(x))
  genes[hasHit] = names(geneRanges)[OLs[xIndex[hasHit],2]]
  
  if ( sum(genes == '?') > 0 & verbose ) catLog('failed', sum(genes == '?'), 'positions...')
  if ( verbose ) catLog('done.\n')

  return(genes)
}


importEnsemblData = function(x, saveDirectory, genome, verbose=T) {
  saveFile = paste0(saveDirectory, '/ensembl', genome, 'annotation.Rdata')
  if ( verbose ) catLog('Linking', length(x), 'positions to genes...')
  chr = xToChr(x, genome)
  pos = xToPos(x, genome)
  
  if ( file.exists(saveFile) ) load(saveFile)
  else {
    if ( verbose ) catLog('waiting for biomaRt/ensembl server...')
    if ( genome == 'hg19' ) {
      mart = useMart(biomart='ensembl', dataset = 'hsapiens_gene_ensembl',
        version='Ensembl Genes', host='grch37.ensembl.org')
      symbolName = 'hgnc_symbol'
      ensemblName = 'ensembl_transcript_id'
      biotype = 'gene_biotype'
    }
    else if ( genome == 'mm10' ) {
      mart = useMart(biomart='ensembl', dataset = 'mmusculus_gene_ensembl')
      symbolName = 'mgi_symbol'
      ensemblName = 'ensembl_transcript_id'
      biotype = 'gene_biotype'
    }
    else stop('hg19 and mm10 are the only supported genomes atm, sorry. :(')
    
    maxLength = 1000
    if ( length(chr) > maxLength ) {
      if ( verbose ) catLog('by batch of', maxLength, '.')
      bm = getBM(attributes=c('chromosome_name', 'start_position', 'end_position', symbolName, ensemblName, biotype), filters = c('chromosome_name', 'start', 'end'),
        value=list(chr[1:maxLength], pos[1:maxLength]-1000, pos[1:maxLength]+1000), mart=mart)
      chr = chr[-(1:maxLength)]
      pos = pos[-(1:maxLength)]

      while(length(chr) > maxLength) {
        if ( verbose ) catLog('.')
        moreBm = getBM(attributes=c('chromosome_name', 'start_position', 'end_position', symbolName, ensemblName, biotype), filters = c('chromosome_name', 'start', 'end'),
          value=list(chr[1:maxLength], pos[1:maxLength]-1000, pos[1:maxLength]+1000), mart=mart)
        chr = chr[-(1:maxLength)]
        pos = pos[-(1:maxLength)]
        bm = rbind(bm, moreBm)
      }

      if ( verbose ) catLog('.')
      moreBm = getBM(attributes=c('chromosome_name', 'start_position', 'end_position', symbolName, ensemblName, biotype), filters = c('chromosome_name', 'start', 'end'),
        value=list(chr, pos-1000, pos+1000), mart=mart)
      bm = rbind(bm, moreBm)
      chr = xToChr(x, genome)
      pos = xToPos(x, genome)
    }
    else {
      bm = getBM(attributes=c('chromosome_name', 'start_position', 'end_position', symbolName, ensemblName, biotype), filters = c('chromosome_name', 'start', 'end'), value=list(chr, pos-1000, pos+1000), mart=mart)
    }
    
    #make sure any empty symbol entries are last, so that we can default to first hit
    bm = bm[order(bm[[symbolName]], decreasing=T),]
    #replace empty symbol entries with ensemble code and remove duplicates
    bm[[symbolName]][bm[[symbolName]] == ''] = bm[[ensemblName]][bm[[symbolName]] == '']
    bm = bm[!duplicated(bm[[symbolName]]),]
    
    if ( file.exists(dirname(saveFile)) ) save(bm, file=saveFile)
  }

  return(bm)
}


hasParalog = function(genes, genome) {
  
  if ( genome == 'hg19' )
    mart = useMart(biomart='ensembl', dataset = 'hsapiens_gene_ensembl',
      version='Ensembl Genes', host='grch37.ensembl.org')
  bm = getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'gene_biotype'), filters = c('hgnc_symbol'), value=list(genes), mart=mart)

  ensemblGenes = unique(bm$ensembl_transcript_id)
  paralogs = getBM(attributes=c('hsapiens_paralog_ensembl_gene', 'ensembl_transcript_id', 'hsapiens_paralog_paralogy_confidence'), filters = c('ensembl_transcript_id'), value=list(ensemblGenes), mart=mart)

  paralogs = paralogs[paralogs$hsapiens_paralog_ensembl_gene != '',]
  ensemblHasParalog = ensemblGenes %in% paralogs$ensembl_transcript_id
  ensemblWithParalog = ensemblGenes[ensemblHasParalog]
    
  bmWithParalog = bm[bm$ensembl_transcript_id %in% ensemblWithParalog,]
  hasParalog = genes %in% bmWithParalog$hgnc_symbol
}
