

#' Run VEP on the provided variants
#'
#' @param variants variants: The variants to be VEPed.
#' @param plotDir character: The plot directory where outputSomatics have been run.
#' @param cpus integer: The number of cpus to be used as most.
#' @param forceRedoVEP Logical: if VEP should be rerun even if saved data already exists.
#'
#' @details This function calls VEP on the output from outputSomaticVariants. For this, VEP needs to be callable by system('vep').
runVEP = function(variants, plotDir, cpus=1, genome='hg19', forceRedoVEP=F) {
  catLog('VEP-ing..')
  dir = paste0(plotDir, '/somatics/')
  for ( name in names(variants$variants) ) {
    infile = paste0(dir, '/', name, '.txt')
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    if ( !file.exists(VEPfile) | forceRedoVEP ) {
      wd = getwd()
      catLog('Moving to ', dir, '\n')
      setwd(dir)
      catLog(name, ': ', sum(variants$variants[[name]]$somaticP > 0 & !is.na(variants$variants[[name]]$somaticP)), ' variants.\n', sep='')
      if ( sum(variants$variants[[name]]$somaticP > 0 & !is.na(variants$variants[[name]]$somaticP)) == 0 ) {
        catLog('Moving back to ', wd, '\n')
        setwd(wd)
        next
      }
      if ( cpus == 1 )
        call = paste0('vep -i ', basename(infile), ' -o ', basename(VEPfile), ' --everything --force_overwrite')
      else
        call = paste0('vep -i ', basename(infile), ' -o ', basename(VEPfile), ' --everything --force_overwrite --fork ', cpus)
      if ( genome == 'mm10' ) call = paste0(call, ' --species mus_musculus')
      catLog(call, '\n')
      call = paste0(call, '\n')
      systemRet = system(call, intern=T)
      if ( !any(grepl('Finished', systemRet)) ) warning('VEP run didnt finish!')
      catLog('Moving back to ', wd, '\n')
      setwd(wd)
    }
  }
  catLog('done.\n')
  
  catLog('Importing VEP results:\n')
  for ( name in names(variants$variants) ) {
    catLog(name, ': ', sep='')
    if ( sum(variants$variants[[name]]$somaticP > 0 & !is.na(variants$variants[[name]]$somaticP)) == 0 ) {
      catLog('no somatic variants.\n')
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
      next
    }
    VEPfile = paste0(dir, '/', name, '.VEP.txt')
    VEPdata = try(read.table(VEPfile, fill=T, header=F), silent=T)
    if ( class(VEPdata) == 'try-error' ) {
      catLog('failed to read VEP file ', VEPfile,'.\n')
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
    }
    else {
      ID = as.numeric(as.factor(VEPdata$V1))
      chr = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[1])
      pos = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[2])
      end = as.numeric(gsub('^.*-', '',pos))
      pos = as.numeric(gsub('-.*$', '',pos))
      var = as.character(VEPdata$V3)
      deletions = grepl('-', var)
      if ( any(deletions) ) pos[deletions] = pos[deletions] - 1
      type = as.character(VEPdata$V7)
      sev = sapply(type, typeToSeverity)
      x = chrToX(chr, pos, genome=genome)

      lastCol = strsplit(as.character(VEPdata$V14), ';')
      polyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) ) return(as.numeric(gsub('\\)$', '', gsub('^.*\\(', '', strs[grep('PolyPhen=',strs)[1]]))))
        else return(0.05)
        })
      
      namevar = var
      namevar[nchar(namevar) > 1] = paste0('+', substring(var[nchar(var) > 1], 2))
      namevar[namevar == '-'] = paste0('-', pmax(1, as.numeric(end)[namevar == '-']-as.numeric(pos)[namevar == '-']))
      rowNames = paste0(x, namevar)
      qNames = rep('', length(unique(ID)))
      qNames[ID] = rowNames
      
      #find most severe effect
      mostSev = rep('unkown', length(unique(ID)))
      sevScore = rep(100, length(unique(ID)))
      for ( i in 1:nrow(VEPdata) ) {
        sevI = sev[i]
        IDI = ID[i]
        if ( sevI - polyPhen[i] < sevScore[IDI] ) {
          sevScore[IDI] = sevI - polyPhen[i]
          mostSev[IDI] = severityToType(sevI)
        }
      }
      
      #add effect and severity to q
      variants$variants[[name]]$type = rep('notChecked', nrow(variants$variants[[name]]))
      variants$variants[[name]]$severity = rep(100, nrow(variants$variants[[name]]))
      
      variants$variants[[name]][qNames,]$type = mostSev
      variants$variants[[name]][qNames,]$severity = sevScore
      catLog(length(mostSev), 'VEPed variants.\n')
    }
  }
  catLog('done.\n')
  return(variants)
}

#' Ranks the variant effects
#'
#' @param type character: The effect.
#'
#' @details Returns a low score for severe effects, higher for less severe, 100 for unknown. 11 and below alters the protein.
typeToSeverity = function(type) {
  #first change some (older? newer?) notation into the ensemble SO terms
  type = gsub('non_coding_exon_variant', 'non_coding_transcript_exon_variant', type)
  type = gsub('nc_transcript_variant', 'non_coding_transcript_exon_variant', type)

  if ( grepl('transcript_ablation', type) ) return(1)
  if ( grepl('splice_acceptor_variant', type) ) return(2)
  if ( grepl('splice_donor_variant', type) ) return(3)
  if ( grepl('stop_gained', type) ) return(4)
  if ( grepl('frameshift_variant', type) ) return(5)
  if ( grepl('stop_lost', type) ) return(6)
  if ( grepl('initiator_codon_variant', type) ) return(7)
  if ( grepl('transcript_amplification', type) ) return(8)
  if ( grepl('inframe_insertion', type) ) return(9)
  if ( grepl('inframe_deletion', type) ) return(10)
  if ( grepl('missense_variant', type) ) return(11)
  if ( grepl('splice_region_variant', type) ) return(12)
  if ( grepl('incomplete_terminal_codon_variant', type) ) return(13)
  if ( grepl('stop_retained_variant', type) ) return(14)
  if ( grepl('synonymous_variant', type) ) return(15)
  if ( grepl('coding_sequence_variant', type) ) return(16)
  if ( grepl('mature_miRNA_variant', type) ) return(17)
  if ( grepl('5_prime_UTR_variant', type) ) return(18)
  if ( grepl('3_prime_UTR_variant', type) ) return(19)
  if ( grepl('non_coding_transcript_exon_variant', type) ) return(20)
  if ( grepl('intron_variant', type) ) return(21)
  if ( grepl('NMD_transcript_variant', type) ) return(22)
  if ( grepl('non_coding_transcript_variant', type) ) return(23)
  if ( grepl('upstream_gene_variant', type) ) return(24)
  if ( grepl('downstream_gene_variant', type) ) return(25)
  if ( grepl('TFBS_ablation', type) ) return(26)
  if ( grepl('TFBS_amplification', type) ) return(27)
  if ( grepl('TF_binding_site_variant', type) ) return(28)
  if ( grepl('regulatory_region_ablation', type) ) return(29)
  if ( grepl('regulatory_region_amplification', type) ) return(30)
  if ( grepl('regulatory_region_variant', type) ) return(31)
  if ( grepl('feature_elongation', type) ) return(32)
  if ( grepl('feature_truncation', type) ) return(33)
  if ( grepl('intergenic_variant', type) ) return(34)
  
  if ( grepl('unknown', type) ) return(100)
  catLog('Dont know about the mutation type: ', type, '\n')
  return(100)
}

#' The variant effect as function of severity
#'
#' @param severity integer: The severity rank.
#'
#' @details Translated back from severity rank to the variant effect.
severityToType = function(severity) {
  if ( severity == 1 ) return('transcript_ablation')
  if ( severity == 2 ) return('splice_acceptor_variant')
  if ( severity == 3 ) return('splice_donor_variant')
  if ( severity == 4 ) return('stop_gained')
  if ( severity == 5 ) return('frameshift_variant')
  if ( severity == 6 ) return('stop_lost')
  if ( severity == 7 ) return('initiator_codon_variant')
  if ( severity == 8 ) return('transcript_amplification')
  if ( severity == 9 ) return('inframe_insertion')
  if ( severity == 10 ) return('inframe_deletion')
  if ( severity == 11 ) return('missense_variant')
  if ( severity == 12 ) return('splice_region_variant')
  if ( severity == 13 ) return('incomplete_terminal_codon_variant')
  if ( severity == 14 ) return('stop_retained_variant')
  if ( severity == 15 ) return('synonymous_variant')
  if ( severity == 16 ) return('coding_sequence_variant')
  if ( severity == 17 ) return('mature_miRNA_variant')
  if ( severity == 18 ) return('5_prime_UTR_variant')
  if ( severity == 19 ) return('3_prime_UTR_variant')
  if ( severity == 20 ) return('non_coding_transcript_exon_variant')
  if ( severity == 21 ) return('intron_variant')
  if ( severity == 22 ) return('NMD_transcript_variant')
  if ( severity == 23 ) return('non_coding_transcript_variant')
  if ( severity == 24 ) return('upstream_gene_variant')
  if ( severity == 25 ) return('downstream_gene_variant')
  if ( severity == 26 ) return('TFBS_ablation')
  if ( severity == 27 ) return('TFBS_amplification')
  if ( severity == 28 ) return('TF_binding_site_variant')
  if ( severity == 29 ) return('regulatory_region_ablation')
  if ( severity == 30 ) return('regulatory_region_amplification')
  if ( severity == 31 ) return('regulatory_region_variant')
  if ( severity == 32 ) return('feature_elongation')
  if ( severity == 33 ) return('feature_truncation')
  if ( severity == 34 ) return('intergenic_variant')
  if ( severity == 100 ) return('unknown')
  return('unknown')
}



#' Runs VEP on a R and plot directory that superFreq has been run on.
#'
#' @param outputDirectories named list with Rdirectory and plotDirectory.
#' @param inputFiles named list as in analyse.
#' @param genome character: 'hg19' or 'mm10'
#' @param cpus integer. maximum number of parallel threads.
#' @param cosmicDirectory character: path to COSMIC directory. (optional)
#' @param forceRedo boolean: if TRUE, previously saved data is ignored and overwritten.
#'
#' @export
#' @details Runs VEP on the variants, and saves the variants (retrivable with loadData) with the addition
#'          information from VEP, and COSMIC if provided.
postAnalyseVEP = function(outputDirectories, inputFiles=NA, metaData=NA, genome='hg19', cpus=1, cosmicDirectory='', forceRedo=F) {
  Rdirectory = outputDirectories$Rdirectory
  plotDirectory = outputDirectories$plotDirectory
  data = loadData(Rdirectory)
  if ( !('allVariants' %in% names(data)) ) {
    warning('Cant find a saved allVariants.\n')
    return()
  }
  
  if ( genome == 'hg19' & 'cosmic' %in% names(data$allVariants$variants$variants[[1]]) & !forceRedo ) return()
  if ( genome == 'mm10' & 'isCCGD' %in% names(data$allVariants$variants$variants[[1]]) & !forceRedo ) return()

  backup = paste0(Rdirectory, '/allVariantsPreVEP.Rdata')
  if ( !file.exists(backup) ) {
    catLog('Backing up variants (in case the vep run goes wrong) to', backup, '\n')
    allVariantsPreVEP = data$allVariants
    if ( !is.null(allVariantsPreVEP) )
      save('allVariantsPreVEP', file=backup)
  }

  variants = data$allVariants$variants
  if ( !('allVariants' %in% names(data)) & 'variants' %in% names(data) )
    variants = data$variants
  if ( class(metaData) == 'data.frame' )
    sampleMetaData = metaData
  else if ( class(inputFiles) == 'list' )
    sampleMetaData = importSampleMetaData(inputFiles$metaDataFile)
  else
    error('Neither inputFiles nor metaData was defined.')
  normals = as.logical(gsub('YES', 'T', gsub('NO', 'F', sampleMetaData$NORMAL)))
  names = make.names(sampleMetaData$NAME, unique=T)
  names(normals) = names
  individuals = sampleMetaData$INDIVIDUAL
  timePoints = sampleMetaData$TIMEPOINT
  names(timePoints) = names(individuals) = names
  samplePairs = metaToSamplePairs(names, individuals, normals)
  timeSeries = metaToTimeSeries(names, individuals, normals)

  outputSomaticVariants(variants, genome=genome, plotDirectory=plotDirectory, cpus=cpus, forceRedo=forceRedo)
  variants = runVEP(variants, plotDirectory, cpus=cpus, genome=genome, forceRedoVEP=forceRedo)
  variants = getMoreVEPinfo(variants, plotDirectory, genome=genome, cosmicDirectory=cosmicDirectory)
  allVariants = data$allVariants
  allVariants$variants = variants

  allVariantSaveFile = paste0(Rdirectory, '/allVariants.Rdata')
  catLog('Saving fully annotated variants to', allVariantSaveFile, '...')
  save('allVariants', file=allVariantSaveFile)
  catLog('done.\n')
  storiesSaveFile = paste0(Rdirectory, '/stories.Rdata')
  catLog('Replacing fully annotated variants in', storiesSaveFile, '...')
  load(storiesSaveFile)
  stories$variants = variants
  save('stories', file=storiesSaveFile)
  catLog('done.\n')

  outputSomaticVariants(variants, genome, plotDirectory, cpus=cpus, forceRedo=T)
  makeSNPprogressionPlots(variants, timeSeries=timeSeries, normals = normals, plotDirectory=plotDirectory,
                          genome=genome, forceRedo=T, maxRowCluster=2000)
  makeRiverPlots(data$stories$stories, variants, genome=genome, cpus=cpus, plotDirectory=plotDirectory, forceRedo=T)
  makeScatterPlots(variants, samplePairs, timePoints, plotDirectory, genome=genome, cpus=cpus, forceRedo=T)
  makeCloneScatterPlots(variants, data$stories$stories, samplePairs, individuals, timePoints,
                        plotDirectory, genome=genome, cpus=cpus, forceRedo=T)
  outputNewVariants(variants, samplePairs, genome, plotDirectory, cpus=cpus, forceRedo=T)
}



#' Import more information about the variants
#'
#' @param variants variants: The variants
#' @param plotDirectory character: The plot directory from the analyse call.
#'
#' @details Extract almost all the information from the VEP run, and cross-checks with COSMIC data as well.
#'          getCosmicCounts is called by this function, see details of that for requirements.
getMoreVEPinfo = function(variants, plotDirectory, genome='hg19', cosmicDirectory='') {
  dir = paste0(plotDirectory, '/somatics')
  catLog('Importing more VEP info:\n')
  for ( name in names(variants$variants) ) {
    catLog(name, '.. ')
    if ( sum(variants$variants[[name]]$somaticP > 0 & !is.na(variants$variants[[name]]$somaticP)) == 0 ) {
      catLog('no somatic variants.\n')
      variants$variants[[name]] = addNullAnnotation(variants$variants[[name]], genome=genome)
      next
    }
    else {
      VEPfile = paste0(dir, '/', name, '.VEP.txt')
      VEPdata = try(read.table(VEPfile, fill=T, header=F), silent=T)
      if ( class(VEPdata) == 'try-error' )
        catLog('failed to read VEP file ', VEPfile,'.\n')
      
      ID = as.numeric(as.factor(VEPdata$V1))
      var = as.character(VEPdata$V3)
      chr = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[1])
      pos = sapply(strsplit(as.character(VEPdata$V2), ':'), function(v) v[2])
      end = as.numeric(gsub('^.*-', '',pos))
      pos = as.numeric(gsub('-.*$', '',pos))
      deletions = grepl('-', var)
      if ( any(deletions) ) pos[deletions] = pos[deletions] - 1
      type = as.character(VEPdata$V7)
      sev = sapply(type, typeToSeverity)
      x = chrToX(chr, pos, genome)
      AApos = VEPdata$V10
      AAbefore = gsub('\\/.*$', '', VEPdata$V11)
      AAafter = gsub('^.*\\/', '', VEPdata$V11)

      secondLastCol = strsplit(as.character(VEPdata$V13), ',')
      lastCol = strsplit(as.character(VEPdata$V14), ';')
      polyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) ) return(gsub('^PolyPhen=', '', strs[grep('PolyPhen=',strs)[1]]))
        else return('')
        })
      numericPolyPhen = sapply(lastCol, function(strs) {
        if ( any(grepl('PolyPhen=', strs)) )
          return(as.numeric(gsub('\\)$', '', gsub('^.*\\(', '', strs[grep('PolyPhen=',strs)[1]]))))
        else return(0.0)
        })

      symbol = sapply(lastCol, function(strs) {
        if ( any(grepl('SYMBOL=', strs)) ) return(gsub('^SYMBOL=', '', strs[grep('SYMBOL=',strs)[1]]))
        else return('')
        })
      exon = sapply(lastCol, function(strs) {
        if ( any(grepl('EXON=', strs)) ) return(gsub('^EXON=', '', strs[grep('EXON=',strs)[1]]))
        else return('')
      })
      sift = sapply(lastCol, function(strs) {
        if ( any(grepl('SIFT=', strs)) ) return(gsub('^SIFT=', '', strs[grep('SIFT=',strs)[1]]))
        else return('')
      })
      domain = sapply(lastCol, function(strs) {
        if ( any(grepl('DOMAINS=', strs)) ) return(gsub('\\,.*$', '', gsub('^DOMAINS=', '', strs[grep('DOMAINS=',strs)[1]])))
        else return('')
      })

      if ( genome == 'hg19' ) {
        cosmic = sapply(secondLastCol, function(strs) {
          if ( any(grepl('COSM', strs)) ) return(strs[grep('COSM',strs)[1]])
          else return('')
        })
        cosmicCounts = getCosmicCounts(cosmic, cosmicDirectory=cosmicDirectory)
        isCosmicCensus = cosmicCounts$found
        cosmicCounts = getCosmicCounts(cosmic, cosmicDirectory=cosmicDirectory, onlyCensus=F)
        cosmicVariantDensity = cosmicCounts$variantDensity
        cosmicGeneDensity = cosmicCounts$geneDensity
        
        censusDensity = getCosmicCensusDensity(cosmicDirectory=cosmicDirectory)
        noHitCensus = !isCosmicCensus & symbol %in% names(censusDensity) & sev <= 11
        if ( any(noHitCensus) ) {
          isCosmicCensus[noHitCensus] = T
          cosmicGeneDensity[noHitCensus] = censusDensity[symbol[noHitCensus]]
        }
      }
      if ( genome == 'mm10' & FALSE ) {
        CCGDdata = read.table(paste0(cosmicDirectory, '/CCGD_export.csv'), sep=',', header=T, stringsAsFactors=F)

        censusGenes = unique(CCGDdata$Mouse.Symbol)
        isCCGD = symbol %in% censusGenes

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
        
        #q$CCGDstudies = rep(0, nrow(q))
        #q$CCGDcancerTypes = rep('', nrow(q))
        #q$CCGDcosmic = rep(F, nrow(q))
        #q$CCGDcgc = rep(F, nrow(q))
        #q$CCGDranks = rep('', nrow(q))
        

            
            
      }
      
      namevar = var
      namevar[nchar(namevar) > 1] = paste0('+', substring(var[nchar(var) > 1], 2))
      namevar[namevar == '-'] = paste0('-', pmax(1, as.numeric(end)[namevar == '-']-as.numeric(pos)[namevar == '-']))
      rowNames = paste0(x, namevar)
      qNames = rep('', length(unique(ID)))
      qNames[ID] = rowNames
      
      #find most severe effect
      mostSev = rep('unkown', length(unique(ID)))
      sevScore = rep(100, length(unique(ID)))
      polyPhenRet = rep('', length(unique(ID)))
      exonRet = rep('', length(unique(ID)))
      siftRet = rep('', length(unique(ID)))
      AAposRet = rep('', length(unique(ID)))
      AAbeforeRet = rep('', length(unique(ID)))
      AAafterRet = rep('', length(unique(ID)))
      domainRet = rep('', length(unique(ID)))
      cosmicRet = rep('', length(unique(ID)))
      isCosmicCensusRet = rep(F, length(unique(ID)))
      cosmicVariantDensityRet = rep(0, length(unique(ID)))
      cosmicGeneDensityRet = rep(0, length(unique(ID)))
      for ( i in 1:nrow(VEPdata) ) {
        sevI = sev[i]
        IDI = ID[i]
        if ( sevI - numericPolyPhen[i] < sevScore[IDI] ) {
          sevScore[IDI] = sevI - numericPolyPhen[i]
          mostSev[IDI] = severityToType(sevI)

          polyPhenRet[IDI] = polyPhen[i]
          exonRet[IDI] = exon[i]
          siftRet[IDI] = sift[i]
          AAposRet[IDI] = AApos[i]
          AAbeforeRet[IDI] = AAbefore[i]
          AAafterRet[IDI] = AAafter[i]
          domainRet[IDI] = domain[i]
          if ( genome == 'hg19' ) {
            cosmicRet[IDI] = cosmic[i]
            isCosmicCensusRet[IDI] = isCosmicCensus[i]
            cosmicVariantDensityRet[IDI] = cosmicVariantDensity[i]
            cosmicGeneDensityRet[IDI] = cosmicGeneDensity[i]
          }
        }
      }
      
      #add effect and severity to q
      variants$variants[[name]] = addNullAnnotation(variants$variants[[name]], genome=genome)
      
      variants$variants[[name]][qNames,]$severity = sevScore
      variants$variants[[name]][qNames,]$type = mostSev
      variants$variants[[name]][qNames,]$polyPhen = polyPhenRet
      variants$variants[[name]][qNames,]$sift = siftRet
      variants$variants[[name]][qNames,]$exon = exonRet
      variants$variants[[name]][qNames,]$AApos = AAposRet
      variants$variants[[name]][qNames,]$AAbefore = AAbeforeRet
      variants$variants[[name]][qNames,]$AAafter = AAafterRet
      variants$variants[[name]][qNames,]$domain = domainRet
      if ( genome == 'hg19' ) {
        variants$variants[[name]][qNames,]$cosmic = cosmicRet
        variants$variants[[name]][qNames,]$isCosmicCensus = isCosmicCensusRet
        variants$variants[[name]][qNames,]$cosmicVariantMPM = cosmicVariantDensityRet
        variants$variants[[name]][qNames,]$cosmicGeneMPMPB = cosmicGeneDensityRet
      }
        catLog(length(mostSev), 'VEPed variants.\n')
    }
  }
  catLog('done.\n')
  return(variants)
}

addNullAnnotation = function(q, genome='hg19') {
  q$severity = rep(100, nrow(q))
  q$type = rep('unknown', nrow(q))
  q$polyPhen = rep('', nrow(q))
  q$sift = rep('', nrow(q))
  q$exon = rep('', nrow(q))
  q$AApos = rep('', nrow(q))
  q$AAbefore = rep('', nrow(q))
  q$AAafter = rep('', nrow(q))
  q$domain = rep('', nrow(q))
  if ( genome == 'hg19' ) {
    q$cosmic = rep('', nrow(q))
    q$isCosmicCensus = rep(F, nrow(q))
    q$cosmicVariantMPM = rep(0, nrow(q))
    q$cosmicGeneMPMPB = rep(0, nrow(q))
  }
  if ( genome == 'mm10' ) {
    q$isCCGD = rep(F, nrow(q))
    q$CCGDstudies = rep(0, nrow(q))
    q$CCGDcancerTypes = rep('', nrow(q))
    q$CCGDcosmic = rep(F, nrow(q))
    q$CCGDcgc = rep(F, nrow(q))
    q$CCGDranks = rep('', nrow(q))
  }
  return(q)
}

#' internal function
#'
#' @details Just remembers which columns are added by getMoreVEPinfo
moreVEPnames = function(genome='hg19') {
  if ( genome == 'hg19' ) return(c('polyPhen', 'sift', 'exon', 'AApos', 'AAbefore', 'AAafter', 'domain', 'cosmic', 'isCosmicCensus', 'cosmicVariantMPM', 'cosmicGeneMPMPB'))
  if ( genome == 'mm10' ) return(c('polyPhen', 'sift', 'exon', 'AApos', 'AAbefore', 'AAafter', 'domain', 'isCCGD', 'CCGDstudies', 'CCGDcancerTypes', 'CCGDcosmic', 'CCGDcgc', 'CCGDranks'))
}  



#' Retrieves information about cosmic variant IDs
#'
#' @param cosmic character: The cosmic IDs.
#' @param cosmicDirectory character: The directory where the cosmic files CosmicMutantExportCensus.tsv and CosmicMutantExport.tsv are located.
#' @param onlyCensus logical: Restrict to cosmic census genes.
#'
#' @details For each cosmic ID, return whether it is seen in a census gene, and the mutations per million tumors (MPM) for that variant, as well as the average MPM over the gene (mutations per million tumors per basepair, MPMPB).
getCosmicCounts = function(cosmic, cosmicDirectory='/wehisan/general/academic/grp_leukemia_genomics/data/resources/COSMIC', onlyCensus=T) {
  if ( cosmicDirectory == '' ) {
    warning('COSMIC directory not specified. To access more cosmic data, download CosmicMutantExportCensus.tsv and CosmicMutantExport.tsv from cosmic, put them in a directory (dont change the names of the files please), and provide the directory to postAnalyseVEP, getMoreVEPinfo or getCosmicCounts.')
    return(list(cosmic=cosmic, found=rep(F, length(cosmic)),
                variantDensity=rep(0, length(cosmic)), geneDensity=rep(0, length(cosmic))))
  }
  
  if ( onlyCensus ) countsFile = paste0(cosmicDirectory, '/cosmicCounts.Rdata')
  else countsFile = paste0(cosmicDirectory, '/allCosmicCounts.Rdata')
  if ( file.exists(countsFile) ) load(countsFile)
  else {
    if ( onlyCensus ) cosmicVariantsFile = paste0(cosmicDirectory, '/CosmicMutantExportCensus.tsv')
    else cosmicVariantsFile = paste0(cosmicDirectory, '/CosmicMutantExport.tsv')
    if ( !file.exists(cosmicVariantsFile) ) {
      warning(paste0('couldnt find cosmics file at ', cosmicVariantsFile,'. Will return 0 counts for all variants.'))
    return(list(cosmic=cosmic, found=rep(NA, length(cosmic)),
                variantDensity=rep(NA, length(cosmic)), geneDensity=rep(NA, length(cosmic))))
    }
    cosmicData = read.table(cosmicVariantsFile, fill=T, sep='\t', header=T)
    variantCounts = table(cosmicData$Mutation.ID)
    geneCounts = table(cosmicData$Gene.name)
    geneLengths = aggregate(Gene.CDS.length ~ Gene.name, cosmicData, FUN=max)
    variantGene = aggregate(Gene.name ~ Mutation.ID, cosmicData, FUN=function(genes) genes[1])
    rownames(variantGene) = variantGene[,1]
    nTumors = length(unique(cosmicData$ID_tumour))
    geneDensity = geneCounts/geneLengths[,2]/nTumors*1e6
    variantDensity = variantCounts/nTumors*1e6

    cosmicCounts = list(geneCounts=geneCounts, geneLengths=geneLengths, nTumors=nTumors,
                        geneDensity=geneDensity, variantDensity=variantDensity, variantGene=variantGene)
    save(cosmicCounts, file=countsFile)
  }

  #count times the specific variant and gene is hit.
  found = cosmic %in% names(cosmicCounts$variantDensity)
  variantDensity = rep(0, length(cosmic))
  variantDensity[found] = cosmicCounts$variantDensity[cosmic[found]]
  geneDensity = rep(0, length(cosmic))
  geneDensity[found] = cosmicCounts$geneDensity[cosmicCounts$variantGene[cosmic[found],2]]

  return(list(cosmic=cosmic, found=found, variantDensity=variantDensity, geneDensity=geneDensity))
}


getCosmicCensusDensity = function(cosmicDirectory='/wehisan/general/academic/grp_leukemia_genomics/data/resources/COSMIC') {
  if ( cosmicDirectory == '' ) {
    warning('COSMIC directory not specified. To access more cosmic data, download CosmicMutantExportCensus.tsv and CosmicMutantExport.tsv from cosmic, put them in a directory (dont change the names of the files please), and provide the directory to postAnalyseVEP, getMoreVEPinfo or getCosmicCounts.')
    return(c())
  }
  
  countsFile = paste0(cosmicDirectory, '/cosmicCounts.Rdata')
  if ( file.exists(countsFile) ) load(countsFile)
  else {
    cosmicVariantsFile = paste0(cosmicDirectory, '/CosmicMutantExportCensus.tsv')
    if ( !file.exists(cosmicVariantsFile) ) {
      warning(paste0('couldnt find cosmics file at ', cosmicVariantsFile,'. Will return 0 counts for all variants.'))
      return(list(cosmic=cosmic, found=rep(NA, length(cosmic)),
                  variantDensity=rep(NA, length(cosmic)), geneDensity=NA))
    }
    cosmicData = read.table(cosmicVariantsFile, fill=T, sep='\t', header=T)
    variantCounts = table(cosmicData$Mutation.ID)
    geneCounts = table(cosmicData$Gene.name)
    geneLengths = aggregate(Gene.CDS.length ~ Gene.name, cosmicData, FUN=max)
    variantGene = aggregate(Gene.name ~ Mutation.ID, cosmicData, FUN=function(genes) genes[1])
    rownames(variantGene) = variantGene[,1]
    nTumors = length(unique(cosmicData$ID_tumour))
    geneDensity = geneCounts/geneLengths[,2]/nTumors*1e6
    variantDensity = variantCounts/nTumors*1e6

    cosmicCounts = list(geneCounts=geneCounts, geneLengths=geneLengths, nTumors=nTumors,
                        geneDensity=geneDensity, variantDensity=variantDensity, variantGene=variantGene)
    save(cosmicCounts, file=countsFile)
  }

  #count times the specific variant and gene is hit.
  censusGenes = names(cosmicCounts$geneDensity)
  censusGenes = gsub('\\_.*$', '', censusGenes)
  uniqueCensusGenes = unique(censusGenes)
  censusGeneDensity = sapply(uniqueCensusGenes, function(gene) sum(cosmicCounts$geneDensity[censusGenes == gene]))
    
  return(censusGeneDensity)
}
