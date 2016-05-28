
#' Quality controls the variants in the vcfs, using the bams.
#'
#' @import BiocGenerics
#'
#' @importFrom BiocGenerics strand start start<- end end<-
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges %within%
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom biomaRt useMart getBM
#' @importFrom Rsamtools scanBamHeader ScanBamParam scanBamFlag scanBam
#' @importFrom R.oo charToInt
getVariantsByIndividual = function(metaData, captureRegions, genome, BQoffset, dbDir, Rdirectory, plotDirectory, cpus, forceRedo=F) {
  catLog('Using variants by individual.\n')
  saveFile = paste0(Rdirectory, '/variantsBI.Rdata')
  if ( file.exists(saveFile) & !forceRedo ) {
    catLog('Loading saved variants by individual.\n')
    load(file=saveFile)
    return(variantsBI)
  }
  
  variantsBI = list()

  #iterate over individuals
  for ( individual in unique(metaData$INDIVIDUAL) ) {
    catLog('\nVariants for', individual, '\n')
    samples = metaData$NAME[metaData$INDIVIDUAL == individual]
    positionsSaveFile = paste0(Rdirectory, '/positions.', individual,'.Rdata')
    variantsSaveFile = paste0(Rdirectory, '/variants.', individual,'.Rdata')
    if ( file.exists(variantsSaveFile) & !forceRedo ) {
      catLog('Loading saved variants for', individual,'\n')
      load(file=variantsSaveFile)
      variantsBI = c(variantsBI, variants)
      next
    }
    if ( file.exists(positionsSaveFile) & !forceRedo ) {
      catLog('Loading saved positions for', individual,'\n')
      load(file=positionsSaveFile)
    }
    else {
      #extract positions in any vcf of any sample from this individual
      vcfFiles = metaData[samples,]$VCF
      positions = vcfToPositions(vcfFiles, genome=genome)

      #keep only positions within 300bp of a capture region
      paddedCaptureRegions = captureRegions
      start(paddedCaptureRegions) = start(captureRegions) - 300
      end(paddedCaptureRegions) = end(captureRegions) + 300
      positions = inGene(positions, paddedCaptureRegions, genome=genome)
      inCapture = SNP2GRanges(positions, genome=genome) %within% paddedCaptureRegions
      catLog('Keeping ', sum(inCapture), ' out of ', length(inCapture),
             ' (', round(sum(inCapture)/length(inCapture), 3)*100,
             '%) SNVs that are inside capture regions.\n', sep='')
      positions = positions[inCapture,]
      catLog('saving positions...')
      save(positions, file=positionsSaveFile)
      catLog('done.\n')
    }

    #extract variant information over the identified positions
    bamFiles = metaData[samples,]$BAM
    variants = lapply(bamFiles, function(file) {
      QCsnps(pileups=importQualityScores(positions, file, BQoffset, genome=genome, cpus=cpus)[[1]],
             positions=positions, cpus=cpus)
    })
    names(variants) = samples
    variants = lapply(variants, function(q) q[apply(!is.na(q), 1, any),])
    variants = shareVariants(variants)
    catLog('saving variants...')
    save(variants, file=variantsSaveFile)
    catLog('done.\n')

    variantsBI = c(variantsBI, variants)
  }

  #check db SNP for called variants.
  variantsBI = lapply(variantsBI, function(q) q[!is.na(q$x),])
  variantsBI = matchTodbSNPs(variantsBI, dir=dbDir, genome=genome, cpus=cpus)
  variantsBI = lapply(variantsBI, function(q) q[order(q$x, q$variant),])
    

  #generate SNPs, mainly for backwards compatibility
  q = do.call(rbind, variantsBI)
  q = q[!duplicated(q$x),]
  q = q[order(q$x),]
  inGene = xToGene(q$x, captureRegions=captureRegions, genome=genome)
  names(inGene) = q$x
  SNPs = data.frame(x=q$x, chr=xToChr(q$x, genome), start=xToPos(q$x, genome), end=xToPos(q$x, genome),
                    inGene=inGene, reference=q$reference, variant=q$variant, db=q$db)
  variantsBI = lapply(variantsBI, function(q) {
    q$inGene = inGene[as.character(q$x)]
    return(q)
    })

  variantsBI = list(variants=variantsBI, SNPs=SNPs)

  catLog('Saving variants..')
  save(variantsBI, file=saveFile)
  catLog('done.\n')

  plotFrequencyDiagnostics(variantsBI, plotDirectory)
    
  return(variantsBI)
}

plotFrequencyDiagnostics = function(variants, plotDirectory) {
  diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
  if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
  FreqDirectory = paste0(diagnosticPlotsDirectory, '/frequencyDistribution/')
  catLog('Plotting frequency distributions to ', FreqDirectory,'..', sep='')
  if ( !file.exists(FreqDirectory) ) dir.create(FreqDirectory)
  for ( sample in names(variants$variants) ) {
    catLog(sample, '..', sep='')
    use = variants$variants[[sample]]$cov > 0 & !is.na(variants$variants[[sample]]$cov)
    if ( is.na(sum(use)) ) {
      warning('Found NA coverage in', sample, '...')
      use[is.na(use)] = F
    }
    if ( sum(use) > 0 ) {
      cov = noneg(variants$variants[[sample]]$cov[use] + pmax(-0.4, pmin(0.4, rnorm(sum(use), 0, 0.2))))
      var = pmin(cov, noneg(variants$variants[[sample]]$var[use] + pmax(-0.4, pmin(0.4, rnorm(sum(use), 0, 0.2)))))
      png(paste0(FreqDirectory, sample, '-varcov.png'), height=1000, width=2000, res=144)
      plotColourScatter((var/cov), cov, log='y', xlab='f', ylab='coverage', verbose=F, main=sample)
      cleanUse = use & variants$variants[[sample]]$flag == ''
      if ( any(cleanUse) ) {
        cov = noneg(variants$variants[[sample]]$cov[cleanUse] + pmax(-0.4, pmin(0.4, rnorm(sum(cleanUse), 0, 0.2))))
        var = pmin(cov, noneg(variants$variants[[sample]]$var[cleanUse] + pmax(-0.4, pmin(0.4, rnorm(sum(cleanUse), 0, 0.2)))))
        plotColourScatter((var/cov), cov, log='y', xlab='f', ylab='coverage', verbose=F, main=paste0(sample, ', clean variants'))
      }
      dev.off()
    }
    
    use = variants$variants[[sample]]$var > 0
    if ( any(use) ) {
      pdf(paste0(FreqDirectory, sample, '-hist.pdf'), height=7, width=14)
      hist((variants$variants[[sample]]$var/variants$variants[[sample]]$cov)[use],
           breaks=(0:100)/100, col=mcri('blue'), main='all Variants',
           xlab='variant frequency', ylab='number of variants')
      cleanUse = use & variants$variants[[sample]]$flag == ''
      if ( any(cleanUse) ) {
        hist((variants$variants[[sample]]$var/variants$variants[[sample]]$cov)[cleanUse],
             breaks=(0:100)/100, col=mcri('blue'), main='clean Variants',
             xlab='variant frequency', ylab='number of variants')
      }
      cleanDbUse = cleanUse & variants$variants[[sample]]$db
      if ( any(cleanDbUse) ) {
        hist((variants$variants[[sample]]$var/variants$variants[[sample]]$cov)[cleanDbUse],
             breaks=(0:100)/100, col=mcri('blue'), main='clean dbSNP Variants',
             xlab='variant frequency', ylab='number of variants')
      }
      dev.off()
    }
  }
  catLog('done.\n')
}



#helper function that imports the variants from a vcf file.
#takes a vector of vcf files as input, and outputs a data frame with
#chromosome, start, end, x (genomic coordinate), reference, variant and gene
#for each position. Duplicated positions are removed, and positions that fails to
#align with an ensembl gene is marked with '?'.
vcfToPositions = function(files, genome='hg19') {
  files = unique(files)
  #if more than one file, call each file separately and rbind the outputs.
  if ( length(files) > 1 ) {
    catLog('Found', length(files), 'files.', '\n')
    ret = do.call(rbind, lapply(files, function(file) vcfToPositions(file, genome=genome)))
    ret = ret[!duplicated(ret$x),]
    ret = ret[order(ret$x),]
    return(ret)
  }

  catLog('Reading file ', files, '...', sep='')
  raw = read.table(files, fill=T, skip=0, row.names=NULL, header=F, as.is=T)
  raw = raw[!grepl('^#', raw$V1),]
  catLog('done.\nProcessing data...')
  if ( nrow(raw) == 0 ) return(matrix(1, nrow=0, ncol=22))
  chrs = gsub('MT', 'M', gsub('chr', '', as.character(raw$V1)))
  use = chrs %in% names(chrLengths(genome))
  variant = sapply(strsplit(raw$V5[use], split=','), function(parts) c(parts, '')[1])
  ret = data.frame(
    chr = chrs[use],
    start = as.numeric(as.character(raw$V2[use])),
    end = as.numeric(as.character(raw$V2[use])),
    x = chrToX(chrs[use], as.numeric(as.character(raw$V2[use])), genome=genome),
    reference = substr(as.character(raw$V4[use]), 1, 1),
    variant = variant, stringsAsFactors=F)

  if ( any(ret$variant == ret$reference) )
    ret[ret$variant == ret$reference,]$variant = '-1'
  
  ret = ret[!duplicated(ret$x),]
  ret = ret[order(ret$x),]
  
  rownames(ret) = ret$x
  catLog('done.\n')
  
  catLog('Returning data frame of dimension', dim(ret), '\n')
  return(ret)
}

#helper function converting from chr+bp coordinates into a single coordinate x that runs over all chromosomes.
chrToX = function(chr, bp, genome='hg19') {
  prevChrL = c(0, cumsum(chrLengths(genome)))
  names(prevChrL) = c(names(chrLengths(genome)), 'outside')
  return(prevChrL[gsub('chr', '', chr)] + bp)
}

#Helper function that takes a SNP data frame and returns an GRanges object for the SNPs.
SNP2GRanges = function(SNPs, genome=genome) {
  ir = IRanges(start = SNPs$start, end = SNPs$end, names = rownames(SNPs), width = rep(1, nrow(SNPs)))
  return(GRanges(seqnames = SNPs$chr, ranges = ir, seqlengths = chrLengths(genome)))
}

#hepler function that marks the variants in a SNPs object as db or non db SNPs.
matchTodbSNPs = function(variants, dir='~/data/dbSNP', genome='hg19', cpus=1) {

  variants = lapply(variants, function(q) {
    q$db = rep(F, nrow(q))
    q$dbMAF = rep(NA, nrow(q))
    q$dbValidated = rep(NA, nrow(q))
    return(q)
  })
  
  for (chr in names(chrLengths(genome)) ) {
    chr = gsub('^M$', 'MT', chr)
    RsaveFile = paste0(dir,'/ds_flat_ch', chr, '.Rdata')
    if ( !file.exists(RsaveFile) ) {
      flatFile = paste0(dir,'/ds_flat_ch', chr, '.dbSNP')
      if ( !file.exists(flatFile) ) {
        catLog('Chromosome ', chr, ': no SNP file found at', flatFile,'. Marking all as not dbSNP.\n', sep='')
        if ( chr %in% c('M', 'MT') ) {
          warning(' mitochondrial dbSNP file not found!')
          next
        }
        else stop('dbSNP file not found!')
      }
      db = read.table(flatFile, header = T, fill=T)
      catLog('extracting positions..')
      db = db[db$pos != '?',] #without position, the dbSNP is useless
      catLog('saving positions for future use..')
      save('db', file=RsaveFile)
    }
    else {
      catLog('Chromosome ', chr, ': loading db positions..')
      load(file=RsaveFile)
    }
    
    catLog('matching to variant positions..')
    variants = mclapply(variants, function(q) {
      thisChr = which(xToChr(q$x, genome) == chr)
      if ( length(thisChr) == 0 ) return(q)
      varPos = xToPos(q$x)[thisChr] + ifelse(grepl('[-]', q$variant[thisChr]), 1, 0)
      dbQ = db[db$pos %in% varPos,]

      q$db[thisChr] = varPos %in% dbQ$pos

      dbVal = dbQ[dbQ$validated,]
      q$dbValidated[thisChr] = varPos %in% dbVal$pos

      dbQ = dbQ[order(dbQ$pos, -dbQ$MAF),]
      dbQ = dbQ[!duplicated(dbQ$pos),]
      dbMAF = dbQ$MAF
      names(dbMAF) = dbQ$pos
      q$dbMAF[thisChr][q$db[thisChr]] = dbMAF[as.character(varPos)][q$db[thisChr]]
      return(q)
    }, mc.cores=cpus)
    catLog('done.\n')
  }
  
  return(variants)
}





#helper functions that counts reads in favour of reference and any present variant
#at the given locations for the given bam files.
importQualityScores = function(positions, files, BQoffset, genome='hg19', cpus=1) {
  ret=list()
  chr = as.character(positions$chr)
  for ( file in files ) {
    catLog(as.character(Sys.time()), '\n')
    catLog('Importing', length(chr), 'pileups by chr from', file, '\n')
    if ( cpus > 1 ) {
      listRet = mclapply(unique(chr), function(ch) {
        use = chr == ch
        return(getQuality(file, chr[use], positions$start[use], BQoffset, cpus=1))
      }, mc.cores=cpus)
    }
    else {
      listRet = lapply(unique(chr), function(ch) {
        use = chr == ch
        return(getQuality(file, chr[use], positions$start[use], BQoffset, cpus=1))
      })
    }
    catLog('done!\n')
    ret[[file]] = do.call(c, listRet)
  }
  return(ret)
}
getQuality = function(file, chr, pos, BQoffset, cpus=1) {
  catLog(chr[1], '.. ', sep='')
  if ( any(grepl('^chr1', names(scanBamHeader(file)[[1]]$targets))) )
    chr = paste0('chr', chr)
  which = GRanges(chr, IRanges(pos-3, pos+3))
  chr = gsub('chr', '', chr)
  p1 = ScanBamParam(which=which, what=c('pos', 'seq', 'cigar', 'qual', 'mapq', 'strand'),
    flag=scanBamFlag(isSecondaryAlignment=FALSE))
  index = paste0(file, '.bai')
  if ( !file.exists(index) ) {
    index = gsub('.bam$', '.bai', file)
    if ( !file.exists(index) ) {
      catLog('Could not find an index file for', file, '\nI want either', paste0(file, '.bai'), 'or', index)
      stop(paste('Could not find an index file for', file, '\nI want either', paste0(file, '.bai'), 'or', index))
    }
  }
  regionReads = scanBam(file, index=index, param=p1)
  regionReads = lapply(regionReads, function(reads) {
    keep = !is.na(reads$mapq)
    if ( length(keep) == 0 ) return(reads)
    if ( sum(!keep) > 0 ) {
      reads$pos = reads$pos[keep]
      reads$seq = reads$seq[keep]
      reads$cigar = reads$cigar[keep]
      reads$qual = reads$qual[keep]
      reads$mapq = reads$mapq[keep]
      reads$strand = reads$strand[keep]
    }
    return(reads)
  })
  if ( cpus > 1 )
    pileups = mclapply(1:length(pos), function(i) {
      return(readsToPileup(regionReads[[i]], pos[i], BQoffset=BQoffset))}, mc.cores=cpus)
  else
    pileups = lapply(1:length(pos), function(i) readsToPileup(regionReads[[i]], pos[i], BQoffset=BQoffset))
  catLog('(done ', chr[1], ') ', sep='')
  return(pileups)
}
readsToPileup = function(reads, pos, BQoffset) {
  cigars = reads$cigar
  seqI = unlist(lapply(reads$pos, function(readPos) pos-readPos+1))
  call = rep(NA, length(cigars))
  qual = rep(NA, length(cigars))
  mapq = reads$mapq
  strand = reads$strand
  easy = grepl('^[0-9]*M$', cigars)
  call[easy] = substring(reads$seq[easy], seqI[easy], seqI[easy])
  qual[easy] = charToInt(substring(reads$qual[easy], seqI[easy], seqI[easy]))-BQoffset
  if ( length(cigars) > 0 ) {
    NAcigar = is.na(cigars)
    easy[NAcigar] = T
    call[NAcigar] = ''
    qual[NAcigar] = 0
    mapq[NAcigar] = 0
    strand[NAcigar] = '+'
  }
  if ( any(!easy) ) {
    messyCalls = solveCigars(cigars[!easy], reads$seq[!easy], reads$qual[!easy], seqI[!easy])
    seqI[!easy] = as.numeric(sapply(messyCalls, function(mC) mC[3]))

    #check for short softclipped calls, check if most other well-aligned reads agree, if so, keep, otherwise discard.
    if ( any(grepl('[a-z]', unlist(lapply(messyCalls, function(call) call[1])))) ) {
      softClipped = grep('[a-z]', unlist(lapply(messyCalls, function(call) call[1])))
      keepSoftClipped = rep(F, length(softClipped))
      if ( sum(easy) > 2 ) {
        consensus = substring(reads$seq[easy], seqI[easy]-2, seqI[easy]+2)
        consensus = consensus[nchar(consensus) == 5]
        if ( length(consensus) > 2 ) {
          before = names(sort(table(substring(consensus, 1, 2)),decreasing=T)[1])
          after = names(sort(table(substring(consensus, 4, 5)),decreasing=T)[1])
          if ( nchar(before) == 2 & nchar(after) == 2 ) {
            seqis = as.numeric(unlist(lapply(messyCalls, function(call) call[3])))[softClipped]
            keepSoftClipped =
              substring(reads$seq[!easy][softClipped], seqis-2,seqis-1) %in% c(before, gsub('^.', '^', before), '') &
            substring(reads$seq[!easy][softClipped], seqis+1, seqis+2) %in% c(after, gsub('.$', '$', before), '')
            if ( sum(keepSoftClipped) > 0 ) {
              messyCalls[softClipped][keepSoftClipped] = lapply(messyCalls[softClipped][keepSoftClipped], function(call) {
                call[1] = toupper(call[1])
                return(call)
              })
            }
          }
        }
      }
      if ( sum(!keepSoftClipped) > 0 ) {
        messyCalls[softClipped][!keepSoftClipped] = lapply(messyCalls[softClipped][!keepSoftClipped], function(call) {
          call = c('', '', 0)
          return(call)
        })
      }
    }
    
    call[!easy] = unlist(lapply(messyCalls, function(solution) solution[1]))
    qual[!easy] = unlist(lapply(messyCalls, function(solution) {
      if ( solution[2] == '' ) return(0)
      else return(mean(charToInt(unlist(strsplit(solution[2],split=NULL)))))
      })) - BQoffset
  }

  #look for stuttering in case of deletions/insertions.
  #this is removed for now, due to not really being needed and cpu intensive.
  note = rep('', length(call))
  #if ( length(easy) > 2 ) {
  #  stutterLength = 20
  #  consensus = substring(reads$seq, seqI-stutterLength, seqI+stutterLength)
  #  consensus = consensus[nchar(consensus) == 2*stutterLength+1]
  #  if ( length(consensus) > 2 ) {
  #    before = names(sort(table(substring(consensus, 1, stutterLength)),decreasing=T)[1])
  #    after = names(sort(table(substring(consensus, stutterLength+2, 2*stutterLength+1)),decreasing=T)[1])
  #    for ( repCall in c('A', 'T', 'C', 'G') ) {
  #      if ( do.call(paste0, as.list(rep(repCall, stutterLength))) %in% c(before, after) ) {
  #        note[call == repCall | grepl('[+-]', call)] = 'stutter'
  #      }
  #    }
  #  }
  #}


  ret = data.frame('call'=call, 'qual'=qual, 'mapq'=mapq, 'strand'=strand, 'note'=note)
  if ( any(is.na(as.matrix(ret))) ) {
    catLog('NAs in pileup. Data frame was:\n')
    for ( row in 1:nrow(ret) ) catLog(as.matrix(ret[row,]), '\n')
    catLog('cigars are\n')
    for ( row in 1:nrow(ret) ) catLog(reads$cigar[row], as.character(reads$seq[row]), as.character(reads$qual[row]), '\n')
    catLog('NA pos', pos, '\n')
    catLog('seqI', seqI, '\n')
    stop('NAs in pileup.')
  }
  return(ret[ret$call != '' & !is.na(ret$qual) & !is.na(ret$mapq) & !is.na(ret$strand) & !is.na(ret$call),])
}
solveCigars = function(cigars, seqs, quals, seqIs) {
  cigars = gsub('[X=]', 'M', cigars)
  type = strsplit(gsub('^[0-9]+', '', cigars), split='[0-9]+')
  lengths = strsplit(gsub('[MIDNSHP]+$', '', cigars), split='[MIDNSHP]')

  calls = rep('', length(cigars))
  qs = rep('', length(cigars))
  checked = rep(1, length(cigars))
  unSolved = rep(T, length(cigars))
  for ( i in 1:max(unlist(lapply(type, length))) ) {
    che = checked[unSolved]
    sI = seqIs[unSolved]
    tyIp1 = sapply(type[unSolved], function(t) if ( length(t) < i+1 ) '' else t[i+1])
    tyI = sapply(type[unSolved], function(t) if ( length(t) < i ) '' else t[i])
    se = seqs[unSolved]
    le = lengths[unSolved]
    leI = as.numeric(sapply(lengths[unSolved], function(l) l[i]))
    qu = quals[unSolved]
    
    eOR = tyI == ''  #read ends before position
    if ( any(eOR) ) {
      calls[which(unSolved)[eOR]] = ''
      qs[which(unSolved)[eOR]] = ''
      seqIs[which(unSolved)[eOR]] = 0
    }

    iAP = che == sI+1 & tyI == 'I' #insert at position
    if ( any(iAP) ) {
      calls[which(unSolved)[iAP]] = paste0('+', substring(se[iAP], sI[iAP]+1, sI[iAP]+leI[iAP]))
      qs[which(unSolved)[iAP]] = substring(qu[iAP], sI[iAP]+1, sI[iAP]+leI[iAP])
      seqIs[which(unSolved)[iAP]] = sI[iAP]
    }    

    dAP = che == sI+1 & (tyI == 'D' | tyI == 'N') #deletion at position
    if ( any(dAP) ) {
      calls[which(unSolved)[dAP]] = paste0('-', leI[dAP])
      qs[which(unSolved)[dAP]] = substring(qu[dAP], sI[dAP], sI[dAP])
      seqIs[which(unSolved)[dAP]] = sI[dAP]
    }

    dSP = che > sI & !iAP & !dAP & tyI != 'S'  #deletion skipped over position
    if ( any(dSP) ) {
      calls[which(unSolved)[dSP]] = ''
      qs[which(unSolved)[dSP]] = ''      
      seqIs[which(unSolved)[dSP]] = 0      
    }


    nAP = che + leI - 1 >= sI & che <= sI & tyI == 'M' #normal read until position
    fID = che + leI - 1 == sI & (tyIp1 == 'I' | tyIp1 == 'D' | tyIp1 == 'N')    #following indel
    if ( any(nAP & fID) ) {
      checked[which(unSolved)[nAP & fID]] = checked[which(unSolved)[nAP & fID]] + leI[nAP & fID]
    }
    if ( any(nAP & !fID) ) {
      calls[which(unSolved)[nAP & !fID]] = substring(se[nAP & !fID], sI[nAP & !fID], sI[nAP & !fID])
      qs[which(unSolved)[nAP & !fID]] = substring(qu[nAP & !fID], sI[nAP & !fID], sI[nAP & !fID])
      seqIs[which(unSolved)[nAP & !fID]] = sI[nAP & !fID]
    }
    
    sAPs =  che > sI & tyI == 'S' #soft clip over position at start of read
    lt4 = leI < 4 #3 or fewer bp are soft clipped
    if ( any(sAPs & lt4) ) {
      calls[which(unSolved)[sAPs & lt4]] = tolower(substring(se[sAPs & lt4], (sI+leI)[sAPs & lt4], (sI+leI)[sAPs & lt4]))
      qs[which(unSolved)[sAPs & lt4]] = substring(qu[sAPs & lt4], (sI+leI)[sAPs & lt4], (sI+leI)[sAPs & lt4])
      seqIs[which(unSolved)[sAPs & lt4]] = (sI+leI)[sAPs & lt4]
    }
    if ( any(sAPs & !lt4) ) {
      calls[which(unSolved)[sAPs & !lt4]] = ''
      qs[which(unSolved)[sAPs & !lt4]] = ''
      seqIs[which(unSolved)[sAPs & !lt4]] = 0
    }

    sAPe = che <= sI & che +leI > sI & tyI == 'S' & tyIp1 == '' #soft clip over position at end of read
    lt4 = leI < 4 #3 or fewer bp are soft clipped
    if ( any(sAPe & lt4) ) {
      calls[which(unSolved)[sAPe & lt4]] = tolower(substring(se[sAPe & lt4], sI[sAPe & lt4], sI[sAPe & lt4]))
      qs[which(unSolved)[sAPe & lt4]] = substring(qu[sAPe & lt4], sI[sAPe & lt4], sI[sAPe & lt4])
      seqIs[which(unSolved)[sAPe & lt4]] = sI[sAPe & lt4]
    }
    if ( any(sAPe & !lt4) ) {
      calls[which(unSolved)[sAPe & !lt4]] = ''
      qs[which(unSolved)[sAPe & !lt4]] = ''
      seqIs[which(unSolved)[sAPe & !lt4]] = 0
    }
    
    nBP = !nAP & tyI == 'M'    #normal ending before position
    if ( any(nBP) ) {
      checked[which(unSolved)[nBP]] = checked[which(unSolved)[nBP]] + leI[nBP]
    }

    iBP = !iAP & !sAPs & !sAPe & (tyI == 'I' | tyI == 'S')   #insertion before positions
    if ( any(iBP) ) {
      checked[which(unSolved)[iBP]] = checked[which(unSolved)[iBP]] + leI[iBP]
      seqIs[which(unSolved)[iBP]] = seqIs[which(unSolved)[iBP]] + leI[iBP]
    }

    dBP = !dAP & (tyI == 'D' | tyI == 'N')   #deletion before position
    if ( any(dBP) ) {
      seqIs[which(unSolved)[dBP]] = seqIs[which(unSolved)[dBP]] - leI[dBP]
    }

    unSolved[eOR | iAP | dAP | dSP | (nAP & !fID) | sAPs | sAPe] = F
  }

  solutions = lapply(1:length(cigars), function(i) c(calls[i], qs[i], seqIs[i]))
  return(solutions)
}
QCsnps = function(pileups, positions, cpus=10) {
  references = as.character(positions$reference)
  variants = as.character(positions$variant)
  xs = positions$x
  catLog('QCing', length(pileups), 'positons...')
  catLog(as.character(Sys.time()), '\n')
  ret = mclapply(1:length(pileups), function(i) {
    if ( i/10000 == round(i/10000) ) catLog(i/1000, 'k..', sep='')
    return(QCsnp(pileups[[i]], references[i], xs[i], variant='', defaultVariant=variants[i]))
  }, mc.cores=cpus)
  ret = mergeVariantList(ret)
  ret = ret[ret$reference != ret$variant,]
  rownames(ret) = paste0(ret$x, ret$variant)
  catLog('done.\n')
  catLog('Variants:', nrow(ret), '\n')
  catLog('Unflagged :', sum(ret$flag==''), '\n')
  catLog('Unflagged over 0%  freq :', sum(ret$var > 0 & ret$flag==''), '\n')
  catLog('Unflagged over 20% freq :', sum(ret$var > ret$cov/5 & ret$flag==''), '\n')
  catLog('Median coverage over unflagged variants:', median(ret$cov[ret$flag=='']), '\n')
  catLog('Repeat flags:', length(grep('Rep', ret$flag)), '\n')
  catLog('Mapping quality flags:', length(grep('Mq', ret$flag)), '\n')
  catLog('Base quality flags:', length(grep('Bq', ret$flag)), '\n')
  catLog('Strand bias flags:', length(grep('Sb', ret$flag)), '\n')
  catLog('Single variant read flags:', length(grep('Svr', ret$flag)), '\n')
  #catLog('Single reference read flags:', length(grep('Srr', ret$flag)), '\n')
  catLog('Minor variant flags:', length(grep('Mv', ret$flag)), '\n')
  catLog('Stutter flags:', length(grep('St', ret$flag)), '\n')
  return(ret)
}
mergeVariantList = function(variantList) {
  dat = unlist(variantList)
  nrow = length(variantList)
  ret = data.frame(
    'x'=as.numeric(dat[grep('^x',names(dat))]),
    'reference'=dat[grep('^reference',names(dat))],
    'variant'=dat[grep('^variant',names(dat))],
    'cov'=as.numeric(dat[grep('^cov',names(dat))]),
    'ref'=as.numeric(dat[grep('^ref[0-9]*$',names(dat))]),
    'var'=as.numeric(dat[grep('^var[0-9]*$',names(dat))]),
    'pbq'=as.numeric(dat[grep('^pbq',names(dat))]),
    'pmq'=as.numeric(dat[grep('^pmq',names(dat))]),
    'psr'=as.numeric(dat[grep('^psr',names(dat))]),
    'RIB'=as.numeric(dat[grep('^RIB',names(dat))]),
    'flag'=dat[grep('^flag',names(dat))], stringsAsFactors=F)
  return(ret)
}
QCsnp = function(pileup, reference, x, variant='', defaultVariant='') {
  if ( is.atomic(pileup) ) {
    return(data.frame('x'=x, 'reference'=reference, 'variant'='', 'cov'=0, 'ref'=0, 'var'=0,
                      'pbq'=1, 'pmq'=1, 'psr'=1, 'RIB'=0, 'flag'='', stringsAsFactors=F))
  }

  ref = pileup$call == reference
  if ( variant=='' ) {
    if ( !grepl('-[ATCG]', defaultVariant) )
      variant = unique(c(defaultVariant, as.character(pileup$call[!ref])))
    else
      variant = unique(as.character(pileup$call[!ref]))
    variant = variant[variant != '']
    if ( length(variant) > 1 ) {
      ret = do.call(rbind, lapply(variant, function(var) QCsnp(pileup, reference, x, var)))
      #flag non-zero variants that are not the strictly largest frequency
      minorVariants = ret$var > 0 & (ret$var < max(ret$var) | sum(ret$var == max(ret$var)) > 1)
      ret$flag[minorVariants] = paste0(ret$flag[minorVariants], 'Mv')
      return(ret)
    }
    else if ( length(variant) == 1 ) QCsnp(pileup, reference, x, variant)
    else variant = defaultVariant
  }
  var = pileup$call == variant

  flag = ''
  
  #check for mapQ 0 or 1, if too many, flag as repeat region.
  #then remove those reads.
  if ( any(pileup$mapq < 2) ) {
    if ( sum(pileup$mapq < 2)/nrow(pileup) > 0.5 )
      flag = paste0(flag, 'Rep')
    pileup = pileup[pileup$mapq > 1,]
    ref = pileup$call == reference
    var = pileup$call == variant
  }
  
  #compare base quality scores between variant and reference
  pbq = 1 
  if ( sum(ref) > 0 & sum(var) > 0 ) {
    pbq = if ( sum(ref > 0) ) wilcox.test(pileup$qual[var], pileup$qual[ref], exact=F)$p.value else 1
    if ( pbq < 0.01 & mean(pileup$qual[var]) < 30 & mean(pileup$qual[ref]) - mean(pileup$qual[var]) > 10  ) flag = paste0(flag, 'Bq')
    else if ( mean(pileup$qual[var]) < 20 ) flag = paste0(flag, 'Bq')
    else if ( sum(pileup$qual[var] > 30) < 0.1*sum(var) ) flag = paste0(flag, 'Bq')
    #the wilcox test fails if all the scores are the same, returning NA. handle.
    if ( is.na(pbq) ) pbq=1
  }
  
  #compare mapping quality scores between variant and reference
  pmq = 1 
  if ( sum(ref) > 0 & sum(var) > 0 ) {
    pmq = if ( sum(ref > 0) ) wilcox.test(pileup$mapq[var], pileup$mapq[ref], exact=F)$p.value else 1
    if ( is.na(pmq) ) pmq = 1
    if ( pmq < 0.01 & mean(pileup$mapq[var]) < 30 & mean(pileup$mapq[ref]) - mean(pileup$mapq[var]) > 10 ) flag = paste0(flag, 'Mq')
    else if ( mean(pileup$mapq[var]) < 20 ) flag = paste0(flag, 'Mq')
    else if ( sum(pileup$mapq[var] > 30) < 0.1*sum(var) ) flag = paste0(flag, 'Mq')
    #the wilcox test fails if all the scores are the same, returning NA. Handle.
    if ( is.na(pmq) ) pmq=1
  }
  
  #compare strand ratio between variant and reference
  psr = 1 
  if ( sum(ref) > 0 & sum(var) > 0 ) {
    psr = fisher.test(matrix(
      c(sum(pileup$strand[ref] == '+'), sum(pileup$strand[ref] == '-'),
        sum(pileup$strand[var] == '+'), sum(pileup$strand[var] == '-')), nrow=2))$p.value
    if ( psr < 0.001 ) flag = paste0(flag, 'Sb')
  }

  #compare strand ratio between variant and reference
  if ( any(pileup$note == 'stutter') ) {
    flag = paste0(flag, 'St')
  }

  #probabilities that the base call and mapping is correct.
  bqW = 1-10^(-pileup$qual/10)
  mqW = 1-10^(-pileup$mapq/10)
  qW = bqW*mqW
  RIB = sum(1-qW)/length(qW)
  if ( is.na(RIB) ) RIB = 0

  #count weighted reads.
  refN = round(sum(qW[ref]))
  varN = round(sum(qW[var]))
  cov = round(sum(qW))

  if (varN == 1) flag = paste0(flag, 'Svr')
  #if (refN == 1) flag = paste0(flag, 'Srr')
  
  ret = data.frame('x'=x, 'reference'=reference, 'variant'=variant, 'cov'=cov, 'ref'=refN, 'var'=varN,
    'pbq'=pbq, 'pmq'=pmq, 'psr'=psr, 'RIB'=RIB, 'flag'=flag, stringsAsFactors=F)
  if ( any(is.na(as.matrix(ret))) ) {
    catLog('Na in return value of QCsnp. Return value was:\n')
    for ( row in 1:nrow(ret) ) catLog(as.matrix(ret[row,]), '\n')
    catLog('input pileup:\n')
    for ( row in 1:nrow(pileup) ) catLog(as.matrix(pileup[row,]), '\n')
    catLog('input reference:', reference, '\n')
    catLog('input variant:', variant, '\n')
    catLog('input defaultVariant:', defaultVariant, '\n')
    catLog('input x:', x, '\n')
    stop('NA in return value of QCsnp. Abort.')
  }
  return(ret)
}

#helper function that makes sure all the variants are present in all the data frames.
shareVariants = function(variants) {
  catLog('Adding uncalled variants..')      
  common = Reduce(union, lapply(variants, rownames))
  if ( length(common) == 0 ) return()
  variants = lapply(variants, function(q) {
    new = common[!(common %in% rownames(q))]
    if ( length(new) > 0 ) {
      catLog('Found ', length(new), '..', sep='')      
      x = as.numeric(gsub('[AGNTCagtcn+-].*$', '', new))
      newQ = q[q$x %in% x,]
      newQ = newQ[!duplicated(newQ$x),]
      rownames(newQ) = newQ$x
      newQ = newQ[as.character(x),]
      newQ$var = 0
      newQ$flag = ''
      rownames(newQ) = new
      newQ$variant = gsub('^[0-9]+', '', new)
      q = rbind(q, newQ)
      q = q[common,]
    }
    return(q)
  })
  catLog('done!\n')      
  return(variants)
}

#helper function that marks the SNPs with the gene they are in.
inGene = function(SNPs, genes, noHit = NA, genome='hg19') {
  SNPsGR = SNP2GRanges(SNPs, genome=genome)
  hits = findOverlaps(SNPsGR, genes)
  inGene = rep(noHit, length(SNPsGR))
  inGene[queryHits(hits)] = names(genes)[subjectHits(hits)]
  SNPs$inGene = inGene
  return(SNPs)
}



#Takes the sample variants and normal bam files and capture regions.
#return the variant information for the normals on the positions that the samples are called on.
getNormalVariants = function(variants, bamFiles, names, captureRegions, genome, BQoffset, dbDir, normalRdirectory, Rdirectory, plotDirectory, cpus, forceRedoSNPs=F, forceRedoVariants=F) {
  SNPs = variants$SNPs
  variants = variants$variants
  variantsToCheck = unique(do.call(c, lapply(variants, rownames)))

  saveFile = paste0(Rdirectory, '/normalVariantsBI.Rdata')
  if ( file.exists(saveFile) & !forceRedoVariants ) {
    catLog('Loading saved variants by individual.\n')
    load(file=saveFile)
    return(normalVariantsBI)
  }
  

  #iterate over normals
  normalVariantsBI = list()
  for ( i in 1:length(bamFiles) ) {
    name = names(bamFiles)[i]
    bam = bamFiles[i]
    variantsSaveFile = paste0(Rdirectory, '/normalVariants.', name,'.Rdata')
    if ( file.exists(variantsSaveFile) & !forceRedoVariants ) {
      catLog('Loading saved variants for', name,'\n')
      load(file=variantsSaveFile)
      normalVariantsBI[[name]] = q
      next
    }

    preExistingVariantsFile = paste0(normalRdirectory, '/q', name, '.Rdata')
    if ( file.exists(preExistingVariantsFile) ) {
      #load what variants are already called
      catLog('Loading normal variants from', preExistingVariantsFile, '.\n')
      load(preExistingVariantsFile)
      
      #fill in any missing variants
      missingVariants = variantsToCheck[!(variantsToCheck %in% rownames(q))]
      if ( length(missingVariants) > 0 ) {
        missingX = as.numeric(gsub('[+-ATCGN]', '', missingVariants))
        missingSNPs = SNPs[SNPs$x %in% missingX,]

        catLog('Filling in missing normal variants from ', name,'.\n', sep='')
        qNew =
          QCsnps(pileups=importQualityScores(SNPs, bam, BQoffset, genome=genome, cpus=cpus)[[1]],
                 positions=SNPs, cpus=cpus)
        q = rbind(q, qNew)
        q = q[!duplicated(paste0(q$x, q$variant)),]
        q = q[order(q$x, q$variant),]
        q = q[apply(!is.na(q), 1, any),]
  
        #save the union of the calls for future batches
        catLog('Saving normal variants to', preExistingVariantsFile, '..')
        save(q, file=preExistingVariantsFile)
        catLog('done.\n')
      }
    }
    else {
      #extract variant information over the predetermined positions
      catLog('Normal variants from ', name,'.\n', sep='')
      q = QCsnps(pileups=importQualityScores(SNPs, bam, BQoffset, genome=genome, cpus=cpus)[[1]],
        positions=SNPs, cpus=cpus)
      q = q[apply(!is.na(q), 1, any),]

      #save the union of the calls for future batches
      catLog('Saving normal variants to', preExistingVariantsFile, '..')
      save(q, file=preExistingVariantsFile)
      catLog('done.\n')
    }
    
    
    #filter out any new variants in the normal samples that are not previously seen
    q = q[variantsToCheck[variantsToCheck %in% rownames(q)],]
    normalVariantsBI[[name]] = q
    
    catLog('saving variants...')
    save(q, file=variantsSaveFile)
    catLog('done.\n')
  }
  #if any variants in the cancer samples are not seen in the normals, fill up with reference calls.
  normalVariantsBI = shareVariants(c(normalVariantsBI, variants))[1:length(normalVariantsBI)]

  
  #check db SNP for called variants.
  normalVariantsBI = matchTodbSNPs(normalVariantsBI, dir=dbDir, genome=genome, cpus=cpus)
  normalVariantsBI = lapply(normalVariantsBI, function(q) q[order(q$x, q$variant),])

  normalVariantsBI = list('variants'=normalVariantsBI, 'SNPs'=SNPs)
  
  catLog('Saving variants..')
  save(normalVariantsBI, file=saveFile)
  catLog('done.\n')

  plotFrequencyDiagnostics(normalVariantsBI, plotDirectory)
    
  return(normalVariantsBI)
}

