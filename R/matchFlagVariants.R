
#ensures that all variants are present in both samples and normals.
#flags variants with suspicious behaviour in the normals.
#marks the somatic-looking variants in the samples
#ie non-db variants that are not present in the normals and not flagged as suspicious.
matchFlagVariants = function(variants, normalVariants, individuals, normals, genome, Rdirectory, flaggingVersion='new', cpus=1, byIndividual=F, forceRedoMatchFlag=F) {
  saveFile = paste0(Rdirectory, '/allVariants.Rdata')
  if ( file.exists(saveFile) & !forceRedoMatchFlag ) {
    catLog('Loading final version of combined variants.\n')
    load(file=saveFile)
    catLog('Estimating reference bias.\n')
    setVariantLoss(allVariants$normalVariants$variants)
    return(allVariants)
  }
  variants$variants = lapply(variants$variants, function(q) q[!is.na(q$cov),])
  normalVariants$variants = lapply(normalVariants$variants, function(q) q[!is.na(q$cov),])
  variants = matchVariants(variants, normalVariants)
  normalVariants = matchVariants(normalVariants, variants)

  #Normalise coverage to the number of available reads. Assumes minor variants are noise.
  variants$variants = normaliseCoverage(variants$variants)
  normalVariants$variants = normaliseCoverage(normalVariants$variants)

  #remove boring variants
  present = rowsums(sapply(variants$variants, function(q) q$var > q$cov*0.05)) > 0
  catLog('Keeping ', sum(present), ' out of ', length(present),
         ' (', round(sum(present)/length(!present), 3)*100, '%) SNVs that are present at 5% frequency in at least one sample.\n', sep='')
  variants$variants = lapply(variants$variants, function(q) q[present,])
  normalVariants$variants = lapply(normalVariants$variants, function(q) q[present,])
  variants$SNPs = normalVariants$SNPs = variants$SNPs[variants$SNPs$x %in% variants$variants[[1]]$x,]

  #Use the normals to flag variants that are noisy in the normals
  if ( flaggingVersion == 'new' )
    variants = newFlagFromNormals(variants, normalVariants, genome, cpus=cpus)
  else
    variants = flagFromNormals(variants, normalVariants, genome, cpus=cpus)
  

  #mark somatic variants
  variants = markSomatics(variants, normalVariants, individuals, normals, cpus=cpus)

  if ( byIndividual ) {
    catLog('Trimming uninformative variants by individual...')
    variants = trimVariantsByIndividual(variants, individuals)
    catLog('done.\n')
  }
  
  allVariants = list('variants'=variants, 'normalVariants'=normalVariants)
  catLog('Saving final version of combined variants..')
  save('allVariants', file=saveFile)
  catLog('done.\n')
  return(allVariants)
}



#helper function that ensures that the first variants object have all the variants.
matchVariants = function(vs1, vs2) {
  catLog('Matching', nrow(vs1$variants[[1]]), 'against', nrow(vs2$variants[[1]]), 'variants..')
  vs1$variants = matchInternalQs(vs1$variants)
  vs2$variants = matchInternalQs(vs2$variants)
  vs1$variants = matchQs(vs1$variants, vs2$variants)
  vs1$SNPs = shareSNPs(vs1$SNPs, vs2$SNPs)
  vs1$SNPs = vs1$SNPs[order(vs1$SNPs$x, vs1$SNPs$variant),]
  catLog('to', nrow(vs1$variants[[1]]), 'variants.\n')
  return(vs1)
}
matchQs = function(qs1, qs2) {
  vs1Vars = rownames(qs1[[1]])
  vs2Vars = rownames(qs2[[1]])
  newV1 = setdiff(vs2Vars, vs1Vars)
  is = which(vs2Vars %in% newV1)
  newVars = data.frame(
    x = qs2[[1]]$x[is],
    reference = qs2[[1]]$reference[is],
    variant = qs2[[1]]$variant[is],
    cov = rep(0, length(is)),
    ref = rep(0, length(is)),
    var = rep(0, length(is)),
    pbq = rep(1, length(is)),
    pmq = rep(1, length(is)),
    psr = rep(1, length(is)),
    RIB = rep(1, length(is)),
    flag = rep('', length(is)),
    db = qs2[[1]]$db[is],
    dbValidated = if ( 'dbValidated' %in% names(qs2[[1]]) ) qs2[[1]]$dbValidated[is]
    else rep(NA, length(is)),
    dbMAF = if ( 'dbMAF' %in% names(qs2[[1]]) ) qs2[[1]]$dbMAF[is]
    else rep(NA, length(is)),
    somaticP = rep(0, length(is)),
    germline = rep(FALSE, length(is)),
    type = if ( 'type' %in% names(qs2[[1]]) ) qs2[[1]]$type[is]
    else rep('notChecked', length(is)),
    severity = if ( 'type' %in% names(qs2[[1]]) ) qs2[[1]]$severity[is]
    else rep(100, length(is)),
    row.names = newV1)
  
  #make sure the new variants have the columns of qs2 and qs1
  for ( col in setdiff(union(names(qs1[[1]]), names(qs2[[1]])), names(newVars)) ) {
    newVars[[col]] = rep(NA, nrow(newVars))
  }
  #then make sure qs1 have the columns of the new variants (which will include qs2)
  for ( col in setdiff(names(newVars), names(qs1[[1]])) ) {
    qs1 = lapply(qs1, function(q) {q[[col]] = rep(NA, nrow(q)); return(q)})
  }
  #we are now ready to rbind the data frames
  qs1 = lapply(qs1, function(q) rbind(q, newVars))
  qs1 = lapply(qs1, function(q) q[order(q$x, q$variant),])
  return(qs1)
}

#makes sure all variant data frames in the list qs has all the rows
#if not available, a token entry with coverage 0 is used.
matchInternalQs = function(qs) {
  if ( length(qs) < 2 ) return(qs)
  for ( i in 2:length(qs) ) {
    qs[1] = matchQs(qs[1], qs[i])
  }
  for ( i in 2:length(qs) ) {
    qs[i] = matchQs(qs[i], qs[1])
  }
return(qs)
}


#helper function that flags variants that have suspicious behaviour in the pool of normals.
flagFromNormals = function(variants, normalVariants, genome, cpus=1) {
  #check normals for recurring noise.
  setVariantLoss(normalVariants$variants)
  varN = do.call(cbind, lapply(normalVariants$variants, function(q) q$var))
  covN = do.call(cbind, lapply(normalVariants$variants, function(q) q$cov))
  flags = do.call(cbind, lapply(normalVariants$variants, function(q) q$flag))
  unflagged = rowsums(flags == '') == ncol(flags)
  db = normalVariants$variants[[1]]$db
  f = rowsums(varN)/rowsums(covN)
  f[rowsums(covN) == 0] = 0
  fs = varN/covN
  fs[covN == 0] = 0  
  psN = matrix(pBinom(as.integer(covN), as.integer(varN), rep(f, ncol(covN))), ncol=ncol(covN))
  pSameF = apply(psN, 1, fisherTest)[2,]
  non0 = f > 0.03 & rowsums(varN) > 1
  isLow = fs < 0.1
  isHigh = fs > 0.95
  isHalf = matrix(pBinom(as.integer(covN), as.integer(varN), rep(refBias(0.5), nrow(covN)*ncol(covN))), ncol=ncol(covN)) > 0.01
  consistent = rowsums(isLow | isHigh | isHalf) == ncol(fs)   #are all samples 0, 0.5 or 1?
  normalNoise = (!db & non0 & pSameF > 0.01) | (db & non0 & pSameF > 0.01 & !consistent)
  catLog('Flagged', sum(normalNoise), 'out of', nrow(fs),
         'variants that are recurrently and consistently noisy in normals.\n')

  present = fs > 0.2 & covN >= 10
  variableNormal = (!normalNoise & rowsums(present) > 0 & !db) | (!normalNoise & non0 & !consistent & db)
  catLog('Flagged another', sum(variableNormal), 'out of', nrow(fs),
         'variants that are not db, but significantly present in at least one normal sample.\n')

  meanCov = rowmeans(covN)
  medianCov = median(meanCov[unflagged & meanCov >= 5])
  manyCopies = meanCov > 10*medianCov
  catLog('Flagged another', sum(manyCopies), 'out of', nrow(fs),
         'variants that have a coverage ten times higher than median', medianCov, ' in the normals.\n')
  
  #check if the variants are consistent with the normal noise frequency, and flag Nnc or Nnn.
  variants$variants =
    lapply(variants$variants, function(q) {
        ps = pBinom(q$cov[normalNoise], q$var[normalNoise], f[normalNoise])
        flag = ifelse(ps > 0.01, 'Nnc', 'Nnn')        #normal noise (non)-consistent
        q$flag[normalNoise] = paste0(q$flag[normalNoise], flag)
        vnFlag = ifelse(variableNormal, 'Vn', '')
        q$flag = paste0(q$flag, vnFlag)
        McFlag = ifelse(manyCopies, 'Mc', '')
        q$flag = paste0(q$flag, McFlag)
        return(q)
    })

  #check for non-db SNPs that behave as polymoprhic db SNPs in the normals
  polymorphic = (!db &                                                   #db
                 rowsums(isLow | isHigh | isHalf) == ncol(fs) &          #all consistent
                 rowsums(isLow & !isHalf) > 0 &                          #one strictly ref 
                 rowsums(!isLow & !isHigh & isHalf) > 0 &                #one strictly het
                 pSameF < 0.01)                                          #not same frequency
  
  catLog('Flagged', sum(polymorphic), 'out of', sum(!db),
         ' non-db variants that are consistently polymorphic in normals.\n')

  variants$variants =
    lapply(variants$variants, function(q) {
      ps = pBinom(q$cov[polymorphic], q$var[polymorphic], f[polymorphic])
      q$flag[polymorphic] = paste0(q$flag[polymorphic], 'Pn')  #polymorphic Normal
      return(q)
    })

  catLog('Flagging SNPs that are in noisy regions. New flags by sample:')
  variants$variants =
    lapply(variants$variants, function(q) {
      flagX = q$x[q$flag != '' & q$flag != 'Nr']
      chrs = xToChr(flagX, genome)
      for ( chr in unique(chrs) ) {
        flagXchr = flagX[chrs == chr]
        chrI = which(xToChr(q$x, genome) == chr)
        inNoisyRegion = unlist(mclapply(q$x[chrI], function(x) {
          within1000 = which(abs(flagXchr - x) < 1000)
          within300 = which(abs(flagXchr[within1000] - x) < 300)
          within200 = which(abs(flagXchr[within300] - x) < 200)
          within100 = which(abs(flagXchr[within200] - x) < 100)
          return(length(within100) >= 10 |
                 length(within200) >= 20 |
                 length(within300) >= 30 |
                 length(within1000) >= 50)
        }, mc.cores=cpus))
        q$flag[chrI[inNoisyRegion]] = paste0(q$flag[chrI[inNoisyRegion]], 'Nr')
      }
      catLog(' ', sum(q$flag[chrI] == 'Nr'), sep='')      
      return(q)
    })
  catLog(' done.\n')
  
  return(variants)
}

#new version of flagging from pool of normals. Looks a bit closer and should be able to filter out
#more low frequency crap without taking real variants. Added as option until tested more.
newFlagFromNormals = function(variants, normalVariants, genome, cpus=1) {
  #check normals for recurring noise.
  setVariantLoss(normalVariants$variants)
  varN = do.call(cbind, lapply(normalVariants$variants, function(q) q$var))
  covN = do.call(cbind, lapply(normalVariants$variants, function(q) q$cov))
  flags = do.call(cbind, lapply(normalVariants$variants, function(q) q$flag))
  unflagged = rowsums(flags == '') == ncol(flags)
  db = normalVariants$variants[[1]]$db
  f = rowsums(varN)/rowsums(covN)
  f[rowsums(covN) == 0] = 0
  fs = varN/covN
  fs[covN == 0] = 0  
  psN = matrix(pBinom(as.integer(covN), as.integer(varN), rep(f, ncol(covN))), ncol=ncol(covN))
  pSameF = apply(psN, 1, fisherTest)[2,]
  non0 = f > 0.03 & rowsums(varN) > 1
  isLow = fs < 0.1
  isHigh = fs > 0.95
  isHalf = matrix(pBinom(as.integer(covN), as.integer(varN), rep(refBias(0.5), nrow(covN)*ncol(covN))), ncol=ncol(covN)) > 0.01
  consistent = rowsums(isLow | isHigh | isHalf) == ncol(fs)   #are all samples 0, 0.5 or 1?
  normalNoise = (!db & non0 & pSameF > 0.01) | (db & non0 & pSameF > 0.01 & !consistent)
  catLog('Flagged', sum(normalNoise), 'out of', nrow(fs),
         'variants that are recurrently and consistently noisy in normals.\n')

  present = fs > 0.2 & covN >= 10
  variableNormal = (!normalNoise & rowsums(present) > 0 & !db) | (!normalNoise & rowsums(present) > 0 & !consistent & db)
  catLog('Flagged another', sum(variableNormal), 'out of', nrow(fs),
         'variants that are not db, but significantly present in at least one normal sample.\n')

  meanCov = rowmeans(covN)
  medianCov = median(meanCov[unflagged & meanCov >= 5])
  manyCopies = meanCov > 10*medianCov
  catLog('Flagged another', sum(manyCopies), 'out of', nrow(fs),
         'variants that have a coverage ten times higher than median', medianCov, ' in the normals.\n')

  #variants that are detected in at least one normal, but not enough
  #to be filtered. Filter or not filter these depending on cancer
  #sample frequency.
 detectedMx = (fs > 0.05 | varN > 1) & covN >= 10
 germlineMx = consistent & (isHalf | isHigh)
 noiseDetected = rowsums(detectedMx & !germlineMx) > 0
 noiseBarelyDetected = noiseDetected & !normalNoise

 noiseFs = varN*(!germlineMx)/covN
 noiseFs[covN == 0] = 0
 noiseF = rowsums(varN*(!germlineMx))/rowsums(covN*(!germlineMx))
 noiseF = ifelse(is.na(noiseF), 0, noiseF)

   
  #check if the variants are consistent with the normal noise frequency, and flag Nnc or Nnn.
  catLog('Flagging variants not above normal background level: ')
  variants$variants =
    lapply(variants$variants, function(q) {
        ps = pBinom(q$cov[normalNoise], q$var[normalNoise], f[normalNoise])
        flag = ifelse(ps > 0.01, 'Nnc', 'Nnn')        #normal noise (non)-consistent
        q$flag[normalNoise] = paste0(q$flag[normalNoise], flag)

	qf = q[noiseBarelyDetected,]$var/q[noiseBarelyDetected,]$cov
	isTripleMean = qf > noiseF[noiseBarelyDetected]
	isDoubleMax = qf > apply(noiseFs[noiseBarelyDetected,], 1, max)
	flag = noiseBarelyDetected
	flag[which(noiseBarelyDetected)[isTripleMean & isDoubleMax]] = F
        #Not above background = Nab
        catLog(sum(q$flag=='' & flag), '..', sep='')
        nabFlag = ifelse(flag, 'Nab', '')
        q$flag = paste0(q$flag, nabFlag)
        
        vnFlag = ifelse(variableNormal, 'Vn', '')
        q$flag = paste0(q$flag, vnFlag)
        McFlag = ifelse(manyCopies, 'Mc', '')
        q$flag = paste0(q$flag, McFlag)
        return(q)
    })

  #check for non-db SNPs that behave as polymoprhic db SNPs in the normals
  polymorphic = (!db &                                                   #db
                 rowsums(isLow | isHigh | isHalf) == ncol(fs) &          #all consistent
                 rowsums(isLow & !isHalf) > 0 &                          #one strictly ref 
                 rowsums(!isLow & !isHigh & isHalf) > 0 &                #one strictly het
                 pSameF < 0.01)                                          #not same frequency
  
  catLog('Flagged', sum(polymorphic), 'out of', sum(!db),
         ' non-db variants that are consistently polymorphic in normals.\n')

  variants$variants =
    lapply(variants$variants, function(q) {
      ps = pBinom(q$cov[polymorphic], q$var[polymorphic], f[polymorphic])
      q$flag[polymorphic] = paste0(q$flag[polymorphic], 'Pn')  #polymorphic Normal
      return(q)
    })

  catLog('Flagging SNPs that are in noisy regions. New flags by sample:')
  variants$variants =
    lapply(variants$variants, function(q) {
      flagX = q$x[q$flag != '' & q$flag != 'Nr']
      chrs = xToChr(flagX, genome)
      for ( chr in unique(chrs) ) {
        flagXchr = flagX[chrs == chr]
        chrI = which(xToChr(q$x, genome) == chr)
        inNoisyRegion = unlist(mclapply(q$x[chrI], function(x) {
          within1000 = which(abs(flagXchr - x) < 1000)
          within300 = which(abs(flagXchr[within1000] - x) < 300)
          within200 = which(abs(flagXchr[within300] - x) < 200)
          within100 = which(abs(flagXchr[within200] - x) < 100)
          return(length(within100) >= 10 |
                 length(within200) >= 20 |
                 length(within300) >= 30 |
                 length(within1000) >= 50)
        }, mc.cores=cpus))
        q$flag[chrI[inNoisyRegion]] = paste0(q$flag[chrI[inNoisyRegion]], 'Nr')
      }
      catLog(' ', sum(q$flag[chrI] == 'Nr'), sep='')      
      return(q)
    })
  catLog(' done.\n')
  
  return(variants)
}




#This helper function assign probabilities that variants are somatics.
# the probabilities that the variants are true somatic variants are added in a column 'somaticP'
markSomatics = function(variants, normalVariants, individuals, normals, cpus=cpus) {
  names = names(variants$variants)
  #pair up cancer normals
  correspondingNormal = findCorrespondingNormal(names, individuals, normals)
  CNs = which(!is.na(correspondingNormal))
  names(CNs) = names(variants$variants[CNs])

  #require consistent and very low frequency between normals.
  varN = do.call(cbind, lapply(normalVariants$variants, function(q) q$var))
  covN = do.call(cbind, lapply(normalVariants$variants, function(q) q$cov))
  f = rowsums(varN)/rowsums(covN)
  fs = varN/covN
  fs[covN == 0] = 0  
  f[rowsums(covN) == 0] = 0
  psN = matrix(pBinom(as.integer(covN), as.integer(varN), rep(f, ncol(covN))), ncol=ncol(covN))
  pSameF = apply(psN, 1, fisherTest)[2,]
  
  #count the ratio of non-db SNPs that have the 'Pn' flag, ie polymorphic normal
  observedPolymorphic = length(grep('Pn',variants$variants[[1]]$flag))
  basePairs = 3e9
  nNormals = length(normalVariants$variants)
  polymorphicFrequency = 1-(1-observedPolymorphic/basePairs)^(1/nNormals)
  
  #calculate somatic p-values for the rest of the samples
  for ( name in names ) {
    catLog('Marking somatic mutations in ', name, '..', sep='')
    q = variants$variants[[name]]
    #Require no flag
    use = q$flag == ''
    q = q[use,]
    freq = q$var/q$cov
    freq[is.na(freq)] = -0.02
    normalFreq = f[use]
    
    
    #set p-value from both number of normals (for population-wide frequency uncertainty)
    #and difference in frequency (is the cancer sample really different)
    #should effectively be a cut on how certain it can be from number of normals, even with
    #perfect 0,0,0,0, 0.5 frequencies
    pPolymorphic = 1/(1+nrow(q)/(polymorphicFrequency*basePairs))
    pNormalFreq = pBinom(q$cov, q$var, normalFreq)
    normalOK = pmin(1, noneg((0.05-normalFreq)/0.05))^2*(normalFreq < freq)
    
    if ( !(name %in% names(CNs)) ) catLog('\nWARNING! No matched normal: removing all validated dbSNPs over 0.1% population frequency from the somatic candidates! Remaining somatics will include rare germline variants.\n', sep='')
    commonPopulationVariant = q$db
    if ( 'dbValidated' %in% names(q) ) commonPopulationVariant = commonPopulationVariant & q$dbValidated
    if ( 'dbMAF' %in% names(q) ) commonPopulationVariant = commonPopulationVariant & q$dbMAF > 0.001
    notInNormal = !commonPopulationVariant
    
    if ( name %in% names(CNs) ) {
      catLog('Correcting somatics using', correspondingNormal[name], 'as matched normal.\n')
      qn = variants$variants[[correspondingNormal[name]]][use,]
      referenceNormal = (qn$var <= 0.1*qn$cov)
      pNormalHet = pBinom(qn$cov[referenceNormal], qn$var[referenceNormal], 0.3)
      fdrNormalHet = p.adjust(pNormalHet, method='fdr')
      referenceNormal[referenceNormal] = fdrNormalHet < 0.01
      #check that there is a significant difference in the somatic vs matched normal if anything is seen in the normal
      psameF = unlist(mclapply(which(referenceNormal), function(i)
        fisher.test(matrix(c(q$ref[i], q$var[i], qn$ref[i], qn$var[i]), nrow=2))$p.value,
        mc.cores=cpus))
      referenceNormal[referenceNormal] = psameF < 0.05 & (q$var/q$cov > 0.2 + 2*(qn$var/qn$cov))[referenceNormal]
      notInNormal = referenceNormal
      normalOK = ifelse(qn$cov == 0, 0.5, noneg(1 - 5*qn$var/qn$cov)) #penalty for non-zero normal frequency
    }
    pSampleOK =
      pmax(0.8, p.adjust(q$pbq, method='fdr'))*
        pmax(0.8, p.adjust(q$pmq, method='fdr'))*
          pmax(0.8, p.adjust(q$psr, method='fdr'))
    pZero = p.adjust(dbinom(q$var, q$cov, q$RIB), method='fdr')   #the base quality cut on 30 correpsonds to 0.001 wrong base calls.
    lowFrequencyPenalty = ifelse(q$cov > 0, pmin(1, 5*q$var/q$cov), 0) #penalty below 20% frequency
    lowCoveragePenalty = pmin(1, q$cov/10) #penalty below 10 reads coverage

    somaticP = (1-pPolymorphic)*(1-pNormalFreq)*normalOK*pSampleOK*(1-pZero)*
               notInNormal*lowFrequencyPenalty*lowCoveragePenalty
 
    variants$variants[[name]]$somaticP = 0
    variants$variants[[name]]$somaticP[use] = somaticP
    catLog('got roughly ', sum(variants$variants[[name]]$somaticP > 0.5), ' somatic variants.\n', sep='')
  }

  return(variants)
}

#helper function that matches samples with a normal from the same individual, if present.
findCorrespondingNormal = function(names, individuals, normals) {
  individuals = individuals[names]
  normals = normals[names]
  ret = sapply(names, function(name) {
    ind = individuals[name]
    isCancer = !normals[name]
    hasNormal = any(normals & individuals == ind & names != name)
    if ( !hasNormal ) return(NA)
    return(names[which(normals & individuals == ind & names != name)[1]]) #fix this to support replicate normals!
  })
  names(ret) = names
  return(ret)
}

#calculates the p-value from a two-tailed binomial.
pBinom = function(cov, var, f) {
  if ( length(f) == 1 ) f = rep(f, length(cov))
  use = cov > 0 & var >= 0 & !is.na(cov) & !is.na(var) & f >= 0 & f <= 1 & !is.na(f)
  p = ifelse(!use, NA, 1)
  p[cov==0] = 1
  if ( length(f) != length(cov) | length(var) != length(cov) )
    cat('Length of f', length(f), 'must match length of cov', length(cov), 'or be 1.\n')
  p[use] = pbinom(var[use], cov[use], f[use]) -  dbinom(var[use], cov[use], f[use])/2
  p[use] = 2*pmin(p[use], 1-p[use])
  return(p)
}

#helper function that matches variants between two SNPs objects
shareSNPs = function(SNPs1, SNPs2) {
  SNPs1 = SNPs1[!is.na(SNPs1$x),]
  SNPs1 = SNPs1[!duplicated(SNPs1$x),]
  SNPs2 = SNPs2[!is.na(SNPs2$x),]
  SNPs2 = SNPs2[!duplicated(SNPs2$x),]
  temp = options()$scipen
  options(scipen = 100)
  rownames(SNPs1) = as.character(SNPs1$x)
  rownames(SNPs2) = as.character(SNPs2$x)
  options(scipen = temp)
  newRows = setdiff(rownames(SNPs2), rownames(SNPs1))
  if ( length(newRows) == 0 ) return(SNPs1)
  newSNPs = SNPs2[newRows,colnames(SNPs1)[colnames(SNPs1) %in% colnames(SNPs2)]]
  for ( col in setdiff(names(SNPs2), names(SNPs1)) ) {
    catLog('adding ', col, '..', sep='')
    SNPs1[[col]] = rep(NA, nrow(SNPs1))
  }
  return(rbind(SNPs1, newSNPs))
}

#helper function combining p-values using fishers method.
fisherTest = function(p) {
  p = p[!is.na(p)]
  Xsq = -2*sum(log(p))
  pVal = 1-pchisq(Xsq, df = 2*length(p))
  return(c(Xsq = Xsq, pVal = pVal))
}

setVariantLoss = function(variants, maxLoops = 99, verbose=T) {
  #if called with several samples, take average
  if ( class(variants) == 'list' ) {
    vL = mean(unlist(lapply(variants, function(var) setVariantLoss(var, maxLoops=maxLoops, verbose=F))))
    assign('.variantLoss', vL, envir = .GlobalEnv)
    if ( verbose ) catLog('Average variant loss is', vL, '\n')
    return(.variantLoss)
  }

  rawF = variants$var/variants$cov
  loops = 0
  assign('.variantLoss', 0, envir = .GlobalEnv)
  while(T) {
    loops = loops + 1
    f = refUnbias(rawF)
    use = !is.na(f) & abs(f-0.5) < 0.3 & variants$flag == '' & pBinom(variants$cov, variants$var, refBias(0.5)) > 0.01 & variants$db
    observedMean = sum(variants$var[use])/sum(variants$cov[use])
    if ( sum(use) == 0 ) {
      observedMean = 0.5
      assign('.variantLoss', 0, envir = .GlobalEnv)
      break
    }
    passedVariants = observedMean/(1-observedMean)
    vL = 1 - passedVariants
    if ( vL - variantLoss() < 0.00001 ) {
      if ( verbose ) catLog('Variant loss converged to', vL, 'after', loops, 'iterations.\n')
      assign('.variantLoss', vL, envir = .GlobalEnv)
      break
    }
    assign('.variantLoss', vL, envir = .GlobalEnv)
    if ( loops > maxLoops ) {
      warning('Iterated to find variant loss', loops, 'times without converging. Will not cancel reference bias. This may affect the quality of CNV calls.')
      assign('.variantLoss', 0, envir = .GlobalEnv)
      break
    }
  }
  if ( variantLoss() > 0.25 ) {
    if ( verbose ) catLog('Variant loss estimated to above 25%, which is suspiciously high. Will not correct for reference bias.\n')
    assign('.variantLoss', 0, envir = .GlobalEnv)
  }
  if ( variantLoss() < 0 ) {
    if ( verbose ) catLog('Variant loss estimated to be negative, which is not realistic. Will not correct for reference bias.\n')
    assign('.variantLoss', 0, envir = .GlobalEnv)
  }
  
  return(.variantLoss)
}

#fetches the estimated loss of variants in alignment etc.
variantLoss = function() {
  if ( !exists('.variantLoss') ) {
    assign('.variantLoss', 0, envir = .GlobalEnv)
    warning('.variantLoss not previouly defined. Setting to 0.\n')
  }
  return(get('.variantLoss', envir = .GlobalEnv))
}

#helper functions that handles reference bias.
refBias = function(f, vL=variantLoss()) {
  return(pmin(1, pmax(0, f*(1-vL)/(f*(1-vL) + 1-f))))
}
refUnbias = function(f, vL=variantLoss()) {
  return(pmin(1, pmax(0, f/(1-vL)/(f/(1-vL) + 1-f))))
}
refBiasMirror = function(f, vL=variantLoss()) {
  return(pmin(1, pmax(0, refBias(1 - refUnbias(f, vL), vL))))
}
mirrorDown = function(var, cov, vL=variantLoss()) {
  f = var/cov
  var = ifelse(f > refBias(0.5, vL), cov*refBiasMirror(f, vL), var)
  var[cov == 0] = 0
  return(round(var))
}


normaliseCoverage = function(qs) {
  catLog('Flagging variants with large minor variants: ')
  qs = lapply(qs, function(q) {
    largeMinorVariant = q$cov < 0.8*(q$ref + q$var)
    catLog(sum(largeMinorVariant), '..')
    q$flag[largeMinorVariant] = paste0(q$flag[largeMinorVariant], 'Lmv')
    q$cov = q$ref + q$var
    return(q)
  })
  return(qs)
}


trimVariantsByIndividual = function(variants, individuals) {
  for ( individual in unique(individuals) ) {
    is = names(individuals)[individuals == individual]
    if ( length(is) > 1 )
      checked = apply(sapply(variants$variants[is], function(q) q$cov), 1, sum) > 0
    else
      checked = variants$variants[[is]]$cov > 0
    checked[is.na(checked)] = F
    for ( i in is ) variants$variants[[i]] = variants$variants[[i]][checked,]
  }
  return(variants)
}
