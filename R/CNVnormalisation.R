



callConsistency = function(M, dM, f, df, postHet, call) {

  #estimate clonality from frequency, propagate uncertainty and set weight
  clonalityF = (0.5-f)/(0.5-f+2^callM*(f-callF))
  clonalityErrorF = df*(freqClonality*(2^callM-1)+1)/(0.5-f+2^callM*(f-callF))
  if ( callF == 0.5 ) {
    clonalityF = 0.5
    clonalityErrorF = 1000
  }
  fweight = 1/clonalityErrorF^2

  #estimate clonality from coverage and propagate uncertainty
  clonalityM = (2^M - 1)/(2^callM - 1)
  clonalityErrorM = 2^M*log(2)*dM/abs(2^callM-1)
  if ( callM == 0 ) {
    clonalityM = 0.5
    clonalityErrorM = 1000
  }
  Mweight = 1/clonalityErrorM^2

  #common estimate of cloanlity based on MAF and coverage, propagate uncertainty
  clonality = (clonalityF*fweight + clonalityM*Mweight)/(fweight+Mweight)
  clonality = max(0, min(1, clonality))
  clonalityError = 1/sqrt(1/clonalityErrorF^2 + 1/clonalityErrorM^2)

  #Expected M and f from combined clonality
  fClone = (clonality*(2*callF*2^callM-1) + 1)/(clonality*(2*2^callM - 2) + 2)
  Mclone = log2(1 + (2^callM-1)*clonality)


  sigmaM = abs(M - Mclone)/dM
  sigmaF = abs(f - fClone)/df
  sigma = sqrt((sigmaM^2+sigmaF^2)/2)

  return(c('clonality'=clonality, 'clonalityError'=clonalityError, 'sigma'=sigma))
}


callConsistencyNumeric = function(M, dM, f, df, postHet, call, N = 1000) {
  #add the systematic variance to the setting, to not overestimate confidence in coverage
  dM = sqrt(dM^2 + systematicVariance()^2)

  callF = callTofM(call)['f']
  callM = callTofM(call)['M']
  callM = max(callM, -20)

  clonality = (0:N)/N

  #Expected M and f for the different clonalities
  fClone = (clonality*(2*callF*2^callM-1) + 1)/(clonality*(2*2^callM - 2) + 2)
  Mclone = log2(1 + (2^callM-1)*clonality)

  #calculate how well measured f and M fits with the expected valued
  sigmaM = abs(M - Mclone)/dM
  sigmaF = abs(f - fClone)/df
  sigmaF[fClone == 0.5] = abs(qnorm(postHet/2, 0, 1))
  sigma = sqrt((sigmaM^2+sigmaF^2)/2)

  bestSigma = which(sigma==min(sigma))[1]
  closeEnough = which(sigma <= min(sigma)+1)
  sigma = sigma[bestSigma]
  bestClonality = clonality[bestSigma]
  lowerClonality = clonality[closeEnough[1]]
  higherClonality = clonality[closeEnough[length(closeEnough)]]

  sigma = min(sigma, 3)

  return(c('sigma'=sigma, 'clonality'=bestClonality, 'lowerClonality'=lowerClonality, 'higherClonality'=higherClonality))
}


clusterToConsistencies = function(cluster, shift=0, N=1000) {
  ret = sapply(allCalls(), function(call) {
    f = cluster$f
    ferr = cluster$ferr
    if ( cluster$cov == 0 ) {
      f = 0.5
      ferr = 1000
    }
    callConsistencyNumeric(cluster$M+shift, cluster$width, f, ferr, cluster$postHet, call, N=N)
  })
  return(t(ret))
}

callConsistencies = function(clusters, shift=0, N=1000) {
  ret = lapply(1:nrow(clusters), function(row) {
    clusterToConsistencies(clusters[row,], shift=shift, N=N)
  })
  return(ret)
}


scoreConsistencies = function(consistencies) {
  #convert the prior probabilities of the calls into sigma penalties
  priorSigma = pToSigma(sapply(allCalls(), callPrior))
  priorSigma = priorSigma - min(priorSigma)
  
  #get the sigma of the best call (of any clonality) from each region
  anyClonalityScores = sapply(consistencies, function(consistency) {
    min(priorSigma + consistency[,'sigma'])
    })
    
  #return RMS of the sigmas as score (lower is better fit)
  return(sqrt(mean(anyClonalityScores^2)))
}

pToSigma = function(p) abs(qnorm(p/2, 0, 1))


findShift = function(clusters, nShifts=40, maxShift = 3, cpus=1, plot=T) {
  #set up shifts more densely around 0
  nShiftsHalf = round(nShifts/2)
  shifts = ((0:nShiftsHalf)/nShiftsHalf)^2*maxShift
  shifts = c(-rev(shifts), shifts[-1])

  #get the scores from consistency for each shift
  shiftConsistencies = mclapply(shifts, function(shift) {
    callConsistencies(clusters, shift=shift, N=1000)
  }, mc.cores=cpus)

  #find the shift with best score
  shiftScores = sapply(shiftConsistencies, scoreConsistencies)
  bestShift = shifts[shiftScores == min(shiftScores)][1]

  #return if no shift was best
  if ( bestShift == 0 ) {
    if ( plot ) plotShiftedMAFLFC(clusters)
    return(clusters)
  }

  #if a non-zero shift is better scored, renormalise and rerun to better find local minimum
  #with the denser shifts around 0.
  clusters$M = clusters$M + bestShift
  clusters = findShift(clusters, nShifts=nShifts, maxShift=maxShift, cpus=cpus, plot=plot)
  return(clusters)
}

plotShiftedMAFLFC = function(clusters, shift=0) {
  clusters$M = clusters$M + shift
  return(plotMAFLFC(clusters))
}

