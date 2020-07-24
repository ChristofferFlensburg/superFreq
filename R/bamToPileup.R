


#' Gets pileup from bam over positions.
bamToPileup = function(bam, fasta, positions, index, Rdirectory, BQoffset=33, samtoolsVersion='samtools 1.9') {

  #read in mpileups, split by chromosome if needed.
  chrs = unique(positions$chr)
  #adding a randon number to temp file name in case several processes read from same bam.
  posFile = paste0(Rdirectory, '/', basename(bam), '_pos_temp_', index, '.', sample.int(1e10, 1))
  mpileups = lapply(chrs, function(chr) {
    region = paste0(chr, ':', min(positions$start[positions$chr == chr]), '-',
      max(positions$start[positions$chr == chr]))
    write.table(positions[positions$chr == chr,c('chr', 'start')], file=posFile, sep='\t', quote=F, row.names=F, col.names=F)
    samtoolsCall = paste0('samtools mpileup -OsB -f ', fasta, ' -r ', region, ' -l ', posFile, ' ', bam)
    #cat('\n', samtoolsCall, '\n')
    mpileups = system(samtoolsCall, intern=T, ignore.stderr=T)
    return(mpileups)
  })
  mpileups = do.call(c, mpileups)

  #remove pos file.
  unlink(posFile)

  #check for overflow: lines that don't have the 7 tabs expected
  tabs = sapply(mpileups, function(mp) nchar(mp) - nchar(gsub('\t', '', mp)))
  if ( any(tabs != 7) ) {
    warning('Some positions dont have expected form of the mpileup, removing these. Could be due to very high (~1000+) read depth, or some samtools issues.')
    catLog('\nWARNING: Some positions dont have expected form of the mpileup, removing these. Could be due to very high (~1000+) read depth, or some samtools issues.\n\n')
    mpileups = mpileups[tabs == 7]
  }

  #match positions, and create empty pileups if missing.
  coords = sapply(mpileups, function(mpileup) {
    parts = strsplit(mpileup, split='\t')[[1]]
    return(paste0(parts[1], ':', parts[2]))
  })
  posCoords = paste0(positions$chr, ':', positions$start)
  allMpileups = paste(positions$chr, positions$start, positions$reference, 0, '*', ' ', ' ', sep='\t')
  names(allMpileups) = posCoords
  allMpileups[as.character(coords)] = mpileups

  #parse the mpileup output
  pileups = lapply(allMpileups, function(mpileup) {
    parts = strsplit(mpileup, split='\t')[[1]]

    reference = parts[3]
    string = gsub('((\\^.)|\\$)', '', parts[5])
    base = strsplit(string, split='')[[1]]
    #look for and handle indels
    inss = which(base == '+')
    dels = which(base == '-')
    if ( length(inss) > 0 | length(dels) > 0 ) {
      toRemove = c()
      if ( length(inss) > 0 ) {
        for ( ins in inss ) {
          #insSize = as.numeric(gsub('[a-zA-Z].*$', '', substring(string, ins+1, ins+10)))
          insSize = as.numeric(gsub('[^0-9].*$', '', substring(string, ins+1, ins+10)))
          insertion = substring(string, ins+floor(log10(insSize))+2, ins+floor(log10(insSize))+1+insSize)
          toRemove = c(toRemove, ins-1, (ins+1):(ins+floor(log10(insSize))+1+insSize))
          base[ins] = paste0('+', insertion)
        }
      }
      if ( length(dels) > 0 ) {
        for ( del in dels ) {
          delSize = as.numeric(gsub('[a-zA-Z].*$', '', substring(string, del+1, del+10)))
          deletion = substring(string, del+floor(log10(delSize))+2, del+floor(log10(delSize))+1+delSize)
          toRemove = c(toRemove, del-1, (del+1):(del+floor(log10(delSize))+1+delSize))
          base[del] = paste0('-', delSize, deletion)
        }
      }
      base = base[-toRemove]
    }
	if ( any(grepl('=', base)) ) {
	  warning('Found = in samtools mpileup at ', paste0(parts[[1]], ':', parts[[2]]),'. This can be caused by a samtools bug when reading deletions immediately followed by insertions. Will remove, but this variant will likely not read out properly.')
	  base = gsub('=' , '', base)
	}

    
    strand = ifelse(grepl('[a-z,]', base), '-', '+')
    base[grepl('^-', base)] = gsub('[a-zA-Z]+$', '', base[grepl('^-', base)])
    isRef = base %in% c('.', ',')
    base[isRef] = reference
    base = toupper(base)
    indel = grepl('[+-]', base)
    
    mqColumn = 7
    if ( !is.na(samtoolsVersion) & samtoolsVersion == 'samtools 1.10' ) mqColumn = 8

    baseQuality = charToInt(strsplit(parts[6], split='')[[1]]) - BQoffset
    mapQuality = charToInt(strsplit(parts[mqColumn], split='')[[1]]) - BQoffset

    note = rep('', length(base))
    
    ret = data.frame('call'=base, 'qual'=baseQuality, 'mapq'=mapQuality,
      'strand'=strand, 'note'=note, stringsAsFactors=F)
    ret = ret[ret$call != '*',]
    ret = ret[ret$call != '>',]
    ret = ret[ret$call != '<',]
    return(ret)
    })

  return(pileups)
  #return
}
