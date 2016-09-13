
#Takes a set of vcf files and corresponding bam files and capture regions.
#The function filters variants outside of the capture regions, and a few other filters.
#The function the checks the variants in the bamfiles, and flags if the variant seems suspicious.
#Outputs a data frame for each sample with variant and reference counts, as well as some quality information.
getVariants = function(vcfFiles, bamFiles, names, captureRegions, genome, BQoffset, dbDir, Rdirectory, plotDirectory, cpus, forceRedoSNPs=F, forceRedoVariants=F) {
  SNPsSaveFile = paste0(Rdirectory, '/SNPs.Rdata')
  if ( file.exists(SNPsSaveFile) & !forceRedoSNPs ) {
    catLog('Loading saved SNVs.\n')
    load(file=SNPsSaveFile)
  }
  else {
    catLog('Importing SNV positions..')
    SNPs = try(varScanSNPToMut(vcfFiles, genome=genome))
    if ( class(SNPs) == 'try-error' ) {
      catLog('Error in varScanSNPToMut.\n')
      stop('Error in varScanSNPToMut.')
    }
    if ( any(is.na(as.matrix(SNPs))) ) {
      catLog('NA in SNPs:\n')
      for ( row in which(sapply(1:nrow(SNPs), function(row) any(is.na(SNPs[row,])))) ) catLog(SNPs[row,], '\n')
      stop('NA in SNPs.')
    }
    catLog('done.\n')
    
    catLog('Tagging SNVs within padded capture region..')
    paddedCaptureRegions = captureRegions
    start(paddedCaptureRegions) = start(captureRegions) - 300
    end(paddedCaptureRegions) = end(captureRegions) + 300
    SNPs = inGene(SNPs, paddedCaptureRegions, genome=genome)
    catLog('done.\n')
    
    inCapture = SNP2GRanges(SNPs, genome=genome) %within% paddedCaptureRegions
    catLog('Keeping ', sum(inCapture), ' out of ', length(inCapture),
           ' (', round(sum(inCapture)/length(inCapture), 3)*100, '%) SNVs that are inside capture regions.\n', sep='')
    SNPs = SNPs[inCapture,]

    flag = flagStrandBias(SNPs)
    catLog('Keeping ', sum(!flag), ' out of ', length(flag),
           ' (', round(sum(!flag)/length(flag), 3)*100, '%) SNVs that have consistent strand ratios.\n', sep='')
    SNPs = SNPs[!flag,]
    
    catLog('Matching to dbSNPs.\n')
    SNPs = matchTodbSNPs(SNPs, dir=dbDir, genome=genome)
    
    catLog('Saving SNVs..')
    save(SNPs, file=SNPsSaveFile)
    catLog('done.\n')
  }

  variantsSaveFile = paste0(Rdirectory, '/variants.Rdata')
  if ( file.exists(variantsSaveFile) & !forceRedoVariants ) {
    catLog('Loading saved variants from', variantsSaveFile, '..')
    load(file=variantsSaveFile)
    catLog('done. Loaded variants of dim', dim(variants[[1]]), '\n')
  }
  else {
    catLog('Calculating variants:\n')
    gc()
    variants = lapply(bamFiles, function(file) {
      QCsnps(pileups=importQualityScores(SNPs, file, BQoffset, genome=genome, cpus=cpus)[[1]], SNPs=SNPs, cpus=cpus)})
    names(variants)=names
    variants = lapply(variants, function(q) q[apply(!is.na(q), 1, any),])
    variants = shareVariants(variants)
    
    present = rowSums(sapply(variants, function(q) q$var > q$cov*0.05)) > 0
    catLog('Keeping ', sum(present), ' out of ', length(present),
           ' (', round(sum(present)/length(!present), 3)*100, '%) SNVs that are present at 5% frequency in at least one sample.\n', sep='')
    variants = lapply(variants, function(q) q[present,])
    SNPs = SNPs[SNPs$x %in% variants[[1]]$x,]
    catLog('Saving SNVs..')
    save(SNPs, file=SNPsSaveFile)
    catLog('done.\n')

    for ( i in 1:length(variants) )  variants[[i]]$db = SNPs[as.character(variants[[i]]$x),]$db
    for ( i in 1:length(variants) )  variants[[i]]$dbValidated = SNPs[as.character(variants[[i]]$x),]$dbValidated
    for ( i in 1:length(variants) )  variants[[i]]$dbMAF = SNPs[as.character(variants[[i]]$x),]$dbMAF

    variants = lapply(variants, function(q) q[order(q$x, q$variant),])
    
    catLog('Saving variants..')
    save(variants, file=variantsSaveFile)
    catLog('done.\n')
    
    diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
    if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
    FreqDirectory = paste0(diagnosticPlotsDirectory, '/frequencyDistribution/')
    catLog('Plotting frequency distributions to ', FreqDirectory,'..', sep='')
    if ( !file.exists(FreqDirectory) ) dir.create(FreqDirectory)
    for ( sample in names(variants) ) {
      catLog(sample, '..', sep='')
      use = variants[[sample]]$cov > 0
      cov = variants[[sample]]$cov[use] + 0.2 + noneg(rnorm(sum(use), 0, 0.2)-0.2)
      var = variants[[sample]]$var[use] + 0.2 + noneg(rnorm(sum(use), 0, 0.2)-0.2)
      png(paste0(FreqDirectory, sample, '-varcov.png'), height=10, width=20, res=144, unit='in')
      plotColourScatter(pmin(1,var/cov), cov, log='y', xlab='f', ylab='coverage', verbose=F, main=sample)
      dev.off()
      
      use = variants[[sample]]$var > 0
      if ( any(use) ) {
        pdf(paste0(FreqDirectory, sample, '-hist.pdf'), height=7, width=14)
        hist((variants[[sample]]$var/variants[[sample]]$cov)[use],
             breaks=(0:100)/100, col=mcri('blue'), main='all Variants',
             xlab='variant frequency', ylab='number of variants')
        cleanUse = use & variants[[sample]]$flag == ''
        if ( any(cleanUse) ) {
          hist((variants[[sample]]$var/variants[[sample]]$cov)[cleanUse],
               breaks=(0:100)/100, col=mcri('blue'), main='clean Variants',
               xlab='variant frequency', ylab='number of variants')
        }
        cleanDbUse = cleanUse & variants[[sample]]$db
        if ( any(cleanDbUse) ) {
          hist((variants[[sample]]$var/variants[[sample]]$cov)[cleanDbUse],
               breaks=(0:100)/100, col=mcri('blue'), main='clean dbSNP Variants',
               xlab='variant frequency', ylab='number of variants')
        }
        dev.off()
      }
    }
    catLog('done.\n')  
  }
  
  return(list(SNPs=SNPs, variants=variants))
}








#helper function that imports the variants from a vcf file.
varScanSNPToMut = function(files, genome='hg19') {
  #if more than one file, call each file separately and rbind the outputs.
  if ( length(files) > 1 ) {
    catLog('Found', length(files), 'files.', '\n')
    return(do.call(rbind, lapply(files, function(file) varScanSNPToMut(file, genome=genome))))
  }

  catLog('Reading file ', files, '...', sep='')
  raw = read.table(files, fill=T, skip=0, row.names=NULL, header=F, as.is=T)
  raw = raw[!grepl('^#', raw$V1),]
  catLog('done.\nProcessing data...')
  if ( nrow(raw) == 1 ) return(matrix(1, nrow=0, ncol=22))
  raw = raw[-1,]
  cons = strsplit(as.character(raw$V5), ':')
  strands = strsplit(as.character(raw$V6), ':')
  chrs = gsub('MT', 'M', gsub('chr', '', as.character(raw$V1)))
  ret = data.frame(
    chr = chrs,
    start = as.numeric(as.character(raw$V2)),
    end = as.numeric(as.character(raw$V2)),
    x = chrToX(chrs, as.numeric(as.character(raw$V2)), genome=genome),
    reference = raw$V3,
    variant = raw$V4,
    consensus = as.character(unlist(lapply(cons, function(v) v[1]))),
    reads = as.numeric(unlist(lapply(cons, function(v) v[2]))),
    readsReference = as.numeric(unlist(lapply(cons, function(v) v[3]))),
    readsVariant = as.numeric(unlist(lapply(cons, function(v) v[4]))),
    frequency = as.numeric(gsub('%', '',as.character(unlist(lapply(cons, function(v) v[5])))))/100,
    pValue = as.numeric(unlist(lapply(cons, function(v) v[6]))),
    filter = as.character(unlist(lapply(strands, function(v) v[1]))),
    referencePlus = as.numeric(unlist(lapply(strands, function(v) v[2]))),
    referenceMinus = as.numeric(unlist(lapply(strands, function(v) v[3]))),
    variantPlus = as.numeric(unlist(lapply(strands, function(v) v[4]))),
    variantMinus = as.numeric(unlist(lapply(strands, function(v) v[5]))),
    filterPValue = as.numeric(unlist(lapply(strands, function(v) v[6]))),
    ref = as.numeric(as.character(raw$V7)),
    het = as.numeric(as.character(raw$V8)),
    hom = as.numeric(as.character(raw$V9)),
    nc = as.numeric(as.character(raw$V10)), stringsAsFactors=F)

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

#helper function that flags variants with suspicious strand distribution.
flagStrandBias = function(SNPs) {
  p = sapply(1:nrow(SNPs), function(row) fisher.test(matrix(
    c(SNPs$referenceMinus[row], SNPs$variantMinus[row],
      SNPs$referencePlus[row], SNPs$variantPlus[row]), nrow=2))$p.value)
  p[is.na(p)] = 1
  fdr = p.adjust(p, method='fdr')
  
  flag = fdr < 0.05
  return(flag)
}

#hepler function that marks the variants in a SNPs object as db or non db SNPs.
matchTodbSNPs = function(SNPs, dir='~/data/dbSNP', genome='hg19') {
  SNPs$db = rep(NA, nrow(SNPs))
  SNPs$dbValidated = rep(NA, nrow(SNPs))
  SNPs$dbMAF = rep(NA, nrow(SNPs))
  chrs = names(chrLengths(genome))
  for (chr in chrs ) {
    if ( !any(SNPs$chr == chr) ) {
      catLog('Chromosome ', chr, ': no SNPs found in this chromosome, done.\n', sep='')
      next
    }
    
    RsaveFile = paste0(dir,'/ds_flat_ch', chr, '.Rdata')
    if ( !file.exists(RsaveFile) ) {
      flatFile = paste0(dir,'/ds_flat_ch', chr, '.dbSNP')
      if ( !file.exists(flatFile) ) {
        catLog('Chromosome ', chr, ': no SNP file found at', flatFile,'. Marking all as not dbSNP.\n', sep='')
        stop('dbSNP file not found!')
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
    
    catLog('matching to sample postions..')
    chrSNPsI = which(SNPs$chr == chr)
    varPos = SNPs$start[chrSNPsI] + ifelse(grepl('[-]', SNPs$variant[chrSNPsI]), 1, 0)
    db = db[db$pos %in% varPos,]
    SNPs$db[chrSNPsI] = varPos %in% db$pos

    dbVal = db[db$validated,]
    SNPs$dbValidated[chrSNPsI] = varPos %in% dbVal$pos

    db = db[order(db$pos, -db$MAF),]
    db = db[!duplicated(db$pos),]
    dbMAF = db$MAF
    names(dbMAF) = db$pos
    SNPs$dbMAF[chrSNPsI][SNPs$db[chrSNPsI]] = dbMAF[as.character(varPos)][SNPs$db[chrSNPsI]]
    
    catLog('marked ', sum(SNPs$db[chrSNPsI]), ' positions as dbSNPs.\n')
  }
  return(SNPs)
}





#helper functions that counts reads in favour of reference and any present variant
#at the given locations for the given bam files.
importQualityScores = function(SNPs, files, BQoffset, genome='hg19', cpus=10) {
  ret=list()
  chr = as.character(SNPs$chr)
  for ( file in files ) {
    catLog(as.character(Sys.time()), '\n')
    catLog('Importing', length(chr), 'pileups by chr from', file, '\n')
    if ( cpus > 1 ) {
      listRet = mclapply(unique(chr), function(ch) {
        use = chr == ch
        return(getQuality(file, chr[use], SNPs$start[use], BQoffset, cpus=1))
      }, mc.cores=cpus)
    }
    else {
      listRet = lapply(unique(chr), function(ch) {
        use = chr == ch
        return(getQuality(file, chr[use], SNPs$start[use], BQoffset, cpus=1))
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
  p1 = ScanBamParam(which=which, what=c('pos', 'seq', 'cigar', 'qual', 'mapq', 'strand'),
    flag=scanBamFlag(isSecondaryAlignment=FALSE))
  index = paste0(file, '.bai')
  if ( !file.exists(index) ) {
    index = gsub('.bam$', '.bai', file)
    if ( !file.exists(index) )
      stop(paste('Could not find an index file for', file, '\nI want either', paste0(file, '.bai'), 'or', index))
  }
  regionReads = scanBam(file, index=index, param=p1)
  before = sum(sapply(regionReads, function(reads) length(reads$mapq)))
  regionReads = lapply(regionReads, function(reads) {
    keep = reads$mapq > 0 & !is.na(reads$mapq)
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
  after = sum(sapply(regionReads, function(reads) length(reads$mapq)))
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
  note = rep('', length(call))
  if ( length(easy) > 2 ) {
    stutterLength = 20
    consensus = substring(reads$seq, seqI-stutterLength, seqI+stutterLength)
    consensus = consensus[nchar(consensus) == 2*stutterLength+1]
    if ( length(consensus) > 2 ) {
      before = names(sort(table(substring(consensus, 1, stutterLength)),decreasing=T)[1])
      after = names(sort(table(substring(consensus, stutterLength+2, 2*stutterLength+1)),decreasing=T)[1])
      for ( repCall in c('A', 'T', 'C', 'G') ) {
        if ( do.call(paste0, as.list(rep(repCall, stutterLength))) %in% c(before, after) ) {
          note[call == repCall | grepl('[+-]', call)] = 'stutter'
        }
      }
    }
  }


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
QCsnps = function(pileups, SNPs, cpus=10) {
  references = as.character(SNPs$reference)
  variants = as.character(SNPs$variant)
  xs = SNPs$x
  catLog('QCing', length(pileups), 'positons...')
  catLog(as.character(Sys.time()), '\n')
  ret = mclapply(1:length(pileups), function(i) {
    if ( i/10000 == round(i/10000) ) catLog(i/1000, 'k..', sep='')
    return(QCsnp(pileups[[i]], references[i], xs[i], variant='', defaultVariant=variants[i]))
  }, mc.cores=cpus)
  ret = mergeVariantList(ret)
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
  
  #check for mapQ 0 or 1, if present, flag as repeat region.
  #then remove those reads.
  if ( any(pileup$mapq < 2) ) {
    if ( sum(pileup$mapq < 2)/nrow(pileup) > 0.1 )
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
      if ( 'type' %in% names(newQ) ) newQ$type = 'notChecked'
      if ( 'severity' %in% names(newQ) ) newQ$severity = 100
      if ( 'somaticP' %in% names(newQ) ) newQ$somaticP = 0
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

  variantsSaveFile = paste0(Rdirectory, '/normalVariants.Rdata')
  normalVariantsSaveFile = paste0(normalRdirectory, '/normalVariants.Rdata')
  if ( file.exists(variantsSaveFile) & !forceRedoVariants ) {
    catLog('Loading saved normal variants from', variantsSaveFile, '..')
    varName = load(file=variantsSaveFile)
    if ( varName == 'variants' ) normalVariants = variants
    catLog('done. Loaded variants of dim', dim(normalVariants[[1]]), '\n')
  }
  else {
    if ( file.exists(normalVariantsSaveFile) ) {
      catLog('Loading saved normal variants from', normalVariantsSaveFile, '..')
      varName = load(file=normalVariantsSaveFile)
      if ( varName == 'variants' ) stop('Normal variants saved as @variants. Should be @normalVariants. Delete save file and rerun.')
      catLog('done. Loaded variants of dim', dim(normalVariants[[1]]), '\n')
    }
    else {
      catLog('Calculating variants:\n')
      gc()
      normalVariants = lapply(bamFiles, function(file) {
        QCsnps(pileups=importQualityScores(SNPs, file, BQoffset, genome=genome, cpus=cpus)[[1]], SNPs=SNPs, cpus=cpus)})
      names(normalVariants)=names
      normalVariants = lapply(normalVariants, function(q) q[apply(!is.na(q), 1, any),])
      catLog('Saving normal variants to normal folder..', sep='')
      normalVariants = shareVariants(normalVariants)
      for ( i in 1:length(normalVariants) )
        normalVariants[[i]]$db = SNPs[as.character(normalVariants[[i]]$x),]$db
      for ( i in 1:length(normalVariants) )
        normalVariants[[i]]$dbValidated = SNPs[as.character(normalVariants[[i]]$x),]$dbValidated
      for ( i in 1:length(normalVariants) )
        normalVariants[[i]]$dbMAF = SNPs[as.character(normalVariants[[i]]$x),]$dbMAF
      save(normalVariants, file=normalVariantsSaveFile)
      catLog('done.\n')
    }
    
    catLog('Adding missing variants..')
    normalX = normalVariants[[1]]$x
    variantX = variants[[1]]$x
    newX = variantX[!(variantX %in% normalX)]
    catLog(length(newX), ' new variants..', sep='')
    if ( length(newX) > 0 ) {
      newSNPs = SNPs[SNPs$x %in% newX,]
      newNormalVariants = lapply(bamFiles, function(file) {
        QCsnps(pileups=importQualityScores(newSNPs, file, BQoffset, genome=genome, cpus=cpus)[[1]], SNPs=newSNPs, cpus=cpus)})
      newNormalVariants = shareVariants(newNormalVariants)
      names(newNormalVariants) = names
      newNormalVariants = lapply(newNormalVariants, function(q) q[apply(!is.na(q), 1, any),])
      for ( i in 1:length(newNormalVariants) )
        newNormalVariants[[i]]$db = SNPs[as.character(newNormalVariants[[i]]$x),]$db
      for ( i in 1:length(newNormalVariants) )
        newNormalVariants[[i]]$dbValidated = SNPs[as.character(newNormalVariants[[i]]$x),]$dbValidated
      for ( i in 1:length(newNormalVariants) )
        newNormalVariants[[i]]$dbMAF = SNPs[as.character(newNormalVariants[[i]]$x),]$dbMAF
      missingColumns = colnames(normalVariants[[1]])[!(colnames(normalVariants[[1]]) %in% colnames(newNormalVariants[[1]]))]
      if ( length(missingColumns) > 0 ) {
        NAcolumns = do.call(cbind, lapply(missingColumns, function(colname) data.frame(rep(NA, nrow(newNormalVariants[[1]])))))
        names(NAcolumns) = missingColumns
        newNormalVariants = lapply(newNormalVariants, function(q) cbind(q, NAcolumns))
      }
      normalVariants = lapply(names, function(name) {
        q = rbind(normalVariants[[name]], newNormalVariants[[name]])
        q = q[order(q$x, q$variant),]
      })
      catLog('saving normal variants to normal folder..', sep='')
      save(normalVariants, file=normalVariantsSaveFile)
      catLog('done.\n')
    }
    catLog('Filtering out variants relevant for this batch...')
    present = normalVariants[[1]]$x %in% variants[[1]]$x
    normalVariants = lapply(normalVariants, function(q) q[present,])
    catLog('done.\n')
    
    catLog('Saving variants..')
    save(normalVariants, file=variantsSaveFile)
    catLog('done.\n')
    
    diagnosticPlotsDirectory = paste0(plotDirectory, '/diagnostics')
    if ( !file.exists(diagnosticPlotsDirectory) ) dir.create(diagnosticPlotsDirectory)
    FreqDirectory = paste0(diagnosticPlotsDirectory, '/frequencyDistribution/')
    catLog('Plotting frequency distributions to ', FreqDirectory,'..', sep='')
    if ( !file.exists(FreqDirectory) ) dir.create(FreqDirectory)
    for ( sample in names(normalVariants) ) {
      catLog(sample, '..', sep='')
      use = normalVariants[[sample]]$cov > 0
      cov = normalVariants[[sample]]$cov[use] + 0.2 + noneg(rnorm(sum(use), 0, 0.2)-0.2)
      var = normalVariants[[sample]]$var[use] + 0.2 + noneg(rnorm(sum(use), 0, 0.2)-0.2)
      png(paste0(FreqDirectory, sample, '-varcov.png'), height=10, width=20, res=144, unit='in')
      plotColourScatter(pmin(1,var/cov), cov, log='y', xlab='f', ylab='coverage', verbose=F, main=sample)
      dev.off()
      
      use = normalVariants[[sample]]$var > 0
      if ( any(use) ) {
        pdf(paste0(FreqDirectory, sample, '-hist.pdf'), height=7, width=14)
        hist((normalVariants[[sample]]$var/normalVariants[[sample]]$cov)[use],
             breaks=(0:100)/100, col=mcri('blue'), main='all Variants',
             xlab='variant frequency', ylab='number of variants')
        cleanUse = use & normalVariants[[sample]]$flag == ''
        if ( any(cleanUse) ) {
          hist((normalVariants[[sample]]$var/normalVariants[[sample]]$cov)[cleanUse],
               breaks=(0:100)/100, col=mcri('blue'), main='clean Variants',
               xlab='variant frequency', ylab='number of variants')
        }
        cleanDbUse = cleanUse & normalVariants[[sample]]$db
        if ( any(cleanDbUse) ) {
          hist((normalVariants[[sample]]$var/normalVariants[[sample]]$cov)[cleanDbUse],
               breaks=(0:100)/100, col=mcri('blue'), main='clean dbSNP Variants',
               xlab='variant frequency', ylab='number of variants')
        }
        dev.off()
      }
    }
    catLog('done.\n')
  }
  
  return(list(SNPs=SNPs, variants=normalVariants))
}



