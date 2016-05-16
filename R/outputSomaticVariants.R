

#prints the somatic variants to an excel sheet.
outputSomaticVariants = function(variants, genome, plotDirectory, cpus=cpus, forceRedo=F) {
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
      toReturn = toReturn[order(somaticP[toReturn], decreasing=T)]
      q = q[toReturn,]
      SNPs = variants$SNPs[variants$SNPs$x %in% q$x,]

      start=xToPos(q$x, genome)
      end = xToPos(q$x, genome)
      variant = q$variant
      reference = q$reference
      if ( any(grepl('-', variant)) ) {
        deletions = grepl('-', variant)
        nDel = as.numeric(gsub('-', '', variant[deletions]))
        reference[deletions] = sapply(nDel, function(n) do.call(paste0, as.list(rep('N', n))))
        variant[deletions] = '-'
        start[deletions] = start[deletions]+1
        end[deletions] = end[deletions]+nDel
      }
      if ( any(grepl('\\+', variant)) ) {
        nIns = nchar(variant[grepl('\\+', variant)])-1
        variant[grepl('\\+', variant)] = paste0(substr(reference[grepl('\\+', variant)], 1, 1), gsub('\\+', '', variant[grepl('\\+', variant)]))
      }

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
        dbMAF=if ( 'dbMAF' %in% names(q) ) q$dbMAF else rep('na', nrow(q)),
        dbValidated=if ( 'dbValidated' %in% names(q) ) q$dbValidated else rep('na', nrow(q)),
        row.names=rownames(q))
      if ( all(moreVEPnames(genome=genome) %in% names(q)) ) {
        catLog('adding more VEP info..')
        somatic = cbind(somatic, q[,moreVEPnames(genome=genome)])
      }
      if ( 'severity' %in% names(q) )  ord = order(q$severity + 10*ifelse(is.na(q$germline), 0, q$germline))
      else ord = order(10*q$germline)
      somatic = somatic[ord,]
      somatics[[sample]] = somatic
      if ( nrow(somatic) > 65000 ) warning('Truncating .xls somatic variants: too many rows.')
      XLsomatics[[sample]] = somatic[1:min(nrow(somatic), 65000),]
    }
    names(XLsomatics) = substring(names(XLsomatics), 1, 30)
    catLog('writing to xls...')
    WriteXLS('XLsomatics', outfile)

    catLog('Writing to .csv...')
    for ( sample in names(somatics) ) {
      somatics[[sample]]$sample = rep(sample, nrow(somatics[[sample]]))
    }
    allSomatics = do.call(rbind, somatics)
    write.csv(allSomatics, gsub('.xls$', '.csv', outfile))
    catLog('done!\n')

    vcfDir = paste0(plotDirectory, '/somatics')
    if ( !file.exists(vcfDir) ) dir.create(vcfDir)
    catLog('Outputting to directory ', vcfDir, '..')
    for ( name in names(somatics) ) {
      somatic = somatics[[name]]
      outfile = paste0(vcfDir, '/', name, '.txt')
      catLog(basename(outfile), '..', sep='')
      if ( nrow(somatic) > 0 ) {
        options(scipen=999)
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
