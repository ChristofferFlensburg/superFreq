

#prints the somatic variants to an excel sheet.
outputSomaticVariants = function(variants, genome, plotDirectory, cpus=cpus, forceRedo=F, onlyForVEP=F, rareGermline=T) {
  vcfDir = paste0(plotDirectory, '/somatics')
  if ( !file.exists(vcfDir) ) dir.create(vcfDir)

  outfile = paste0(plotDirectory, '/somaticVariants.xls')
  if ( (!file.exists(outfile) | forceRedo) ) {
    somatics = list()
    XLsomatics = list()
    catLog('Printing somatic variants to ', outfile, '.\n', sep='')
    for ( sample in names(variants$variants) ) {
      catLog(sample, '..', sep='')
      q = variants$variants[[sample]]
      somaticP = q$somaticP
      toReturn = which(somaticP > 0)
      if ( !rareGermline ) toReturn = which(somaticP > 0 & (!q$germline | is.na(q$germline)))
      toReturn = toReturn[order(somaticP[toReturn], decreasing=T)]
      q = q[toReturn,]

      vcfFile = paste0(vcfDir, '/', sample, '.vcf')
      if ( !onlyForVEP ) writeToVCF(q, vcfFile, genome=genome)

      start = as.integer(xToPos(q$x, genome))
      end = as.integer(xToPos(q$x, genome))
      variant = q$variant
      reference = q$reference

      refvarstartendList = convertVariantsToVCF(list(reference, variant, start, end))
      reference = refvarstartendList[[1]]
      variant = refvarstartendList[[2]]
      start = refvarstartendList[[3]]
      end = refvarstartendList[[4]]
      
      somatic = data.frame(
        chr=xToChr(q$x, genome),
        start=start,
        end=end,
        reference=reference,
        variant=variant,
        inGene=gsub('.+:', '', q$inGene),
        severity=if ( 'severity' %in% names(q) ) q$severity else rep('na', nrow(q)),
        effect=if ( 'type' %in% names(q) ) q$type else rep('notChecked', nrow(q)),
        f=q$var/q$cov,
        cov=q$cov,
        ref=q$ref,
        var=q$var,
        flag=q$flag,
        pbq=q$pbq,
        pmq=q$pmq,
        psr=q$psr,
        somaticScore=q$somaticP,
        germlineLike=ifelse(is.na(q$germline), 'na', ifelse(q$germline, 'YES', '')),
        dbSNP=ifelse(q$db, 'dbSNP', ''),
        dbMAF=if ( 'dbMAF' %in% names(q) ) ifelse(q$db, q$dbMAF, '') else rep('na', nrow(q)),
        dbValidated=if ( 'dbValidated' %in% names(q) ) ifelse(q$db, q$dbValidated, '') else rep('na', nrow(q)),
        ExAC= if ( 'exac' %in% names(q) ) ifelse(q$exac, 'ExAC', '') else rep('na', nrow(q)),
        ExAC_AF = if ( 'exacAF' %in% names(q) ) ifelse(q$exac, q$exacAF, '') else rep('na', nrow(q)),
        ExAC_Filter = if ( 'exacFilter' %in% names(q) ) ifelse(q$exac, q$exacFilter, '') else rep('na', nrow(q)),
        row.names=rownames(q), stringsAsFactors=F)
      if ( all(moreVEPnames(genome=genome) %in% names(q)) ) {
        catLog('adding more VEP info..')
        somatic = cbind(somatic, q[,moreVEPnames(genome=genome)])
      }
      if ( all(annotationColumns(genome=genome) %in% names(q)) ) {
        catLog('adding more VariantAnnotation info..')
        somatic = cbind(somatic, q[,annotationColumns(genome=genome)])
      }
      if ( 'severity' %in% names(q) )  ord = order(q$severity + 10*ifelse(is.na(q$germline), 0, q$germline))
      else ord = order(10*q$germline)
      somatic = somatic[ord,]
      somatics[[sample]] = somatic
      if ( nrow(somatic) > 65000 ) warning('Truncating .xls somatic variants: too many rows.')
      XLsomatics[[sample]] = somatic[1:min(nrow(somatic), 65000),]
    }
    names(XLsomatics) = make.unique(substring(names(XLsomatics), 1, 29))
    if ( !onlyForVEP ) catLog('writing to xls...')
    if ( !onlyForVEP ) WriteXLS('XLsomatics', outfile)

    if ( !onlyForVEP ) catLog('Writing to .csv...')
    for ( sample in names(somatics) ) {
      somatics[[sample]]$sample = rep(sample, nrow(somatics[[sample]]))
    }
    allSomatics = do.call(rbind, somatics)
    if ( !onlyForVEP ) write.csv(allSomatics, gsub('.xls$', '.csv', outfile))
    if ( !onlyForVEP ) catLog('done!\n')

    catLog('Outputting to directory ', vcfDir, '..')
    for ( name in names(somatics) ) {
      somatic = somatics[[name]]
      outfile = paste0(vcfDir, '/', name, '.txt')
      catLog(basename(outfile), '..', sep='')
      if ( nrow(somatic) > 0 ) {
        options(scipen=999)

        insertions = grepl('\\+', rownames(somatic))
        if ( any(insertions) ) {
          somatic$reference[insertions] = '-'
          somatic$start[insertions] = somatic$start[insertions]+1
          somatic$variant[insertions] = gsub('\\+', '', somatic$variant[insertions])
        }

        mx = cbind(as.character(somatic$chr),
          as.character(somatic$start), as.character(somatic$end),
          paste0(as.character(somatic$reference), '/', as.character(somatic$variant)),
          rep('+', nrow(somatic)), as.character(1:nrow(somatic)))
        write(t(mx), file=outfile, sep='\t', ncolumns=ncol(mx))
        options(scipen=5)
      }
      else
        write('No somatic mutations detected.', file=outfile, sep='\t')
    }
    catLog('done.\n')
  }
}


#' exports variants to VCF
#'
#' @param q A variant data frame from superFreq.
#' @param vcfFile The path to the output file.
#' @param genome The genome, such as 'hg19', 'hg38' or 'mm10'. Defaults to 'hg19'.
#' @param SNVonly boolean. Set to TRUE to only output SNVs, not indels. Defaults to FALSE.
#' @param addSomaticP boolean. Set to TRUE to include a column with the somaticP score from superFreq. Defaults to FALSE.
#'
#' @details This function outputs superFreq variants to a VCF for access from other software.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data = loadData(Rdirectory)
#' q = data$allVariants$variants$variants$mySample
#' writeToVCF(q, 'mySample.superFreq.vcf')
#' }
writeToVCF = function(q, vcfFile, genome='hg19', snvOnly=F, addSomaticP=F) {
  if ( snvOnly ) q = q[q$variant %in% c('A', 'T', 'C', 'G'),]

  refvarstartendList = convertVariantsToVCF(list(q$reference, q$variant, xToPos(q$x, genome), xToPos(q$x, genome)))
  reference = refvarstartendList[[1]]
  variant = refvarstartendList[[2]]
  start = refvarstartendList[[3]]
  end = refvarstartendList[[4]]
  
  preambula = c('##fileformat=VCFv4.0',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tCOVERAGE\tVARIANTREADS\tVAF')
  if ( addSomaticP )
  preambula = c('##fileformat=VCFv4.0',
    '##reference=', genome,
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
  chrom = xToChr(q$x, genome)
  pos = as.integer(xToPos(q$x, genome))
  ID = rownames(q)
  ref = reference
  alt = variant
  qual = pmin(60, -log10(1-q$somaticP)*10)
  filter = gsub('^$', 'PASS', expandFlags(q$flag))
  info = paste0('DP=', q$cov, ';AF=', ifelse(q$cov > 0, q$var/q$cov, 0))

  out = cbind(chrom, pos, ID, ref, alt, qual, filter, info)
  if ( addSomaticP )
      out = cbind(chrom, pos, ID, ref, alt, qual, filter, info, q$somaticP)
  body = apply(out, 1, function(strs) do.call(paste, c(as.list(strs), sep='\t')))
  body = do.call(paste, c(as.list(preambula), as.list(body), sep='\n'))

  write.table(body, file=vcfFile, row.names=F, col.names=F, quote=F)
}

expandFlags = function(flags) {
  flags = gsub('(Rep)([A-Z]|$)', 'Repeatregion:\\2', flags)
  flags = gsub('(Bq)([A-Z]|$)', 'Basecallquality:\\2', flags)
  flags = gsub('(Mq)([A-Z]|$)', 'Mappingquality:\\2', flags)
  flags = gsub('(Svr)([A-Z]|$)', 'Singlevariantread:\\2', flags)
  flags = gsub('(Nnc)([A-Z]|$)', 'Normalnoiseconsistent:\\2', flags)
  flags = gsub('(Nnn)([A-Z]|$)', 'Normalnoisenonconsistent:\\2', flags)
  flags = gsub('(Vn)([A-Z]|$)', 'Variablenormals:\\2', flags)
  flags = gsub('(Mv)([A-Z]|$)', 'Minorvariant:\\2', flags)
  flags = gsub('(St)([A-Z]|$)', 'Stutter:\\2', flags)
  flags = gsub('(Sb)([A-Z]|$)', 'Strandbias:\\2', flags)
  flags = gsub('(Pa)([A-Z]|$)', 'Palindromic:\\2', flags)
  flags = gsub('(Mc)([A-Z]|$)', 'Manycopies:\\2', flags)
  flags = gsub('(Pn)([A-Z]|$)', 'Polymorphicnormals:\\2', flags)
  flags = gsub('(Nr)([A-Z]|$)', 'Noisyregion:\\2', flags)

  return(flags)
}
