



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
  dM = sqrt(dM^2 + getSystematicVariance()^2)

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
    callScore = superFreq:::callConsistencyNumeric(cluster$M+shift, cluster$width, f, ferr, cluster$postHet, call, N=N)
    #penalty for large regions of complete loss.
    if ( call == 'CL' & cluster$x2-cluster$x1 > 5e6 ) callScore['sigma'] = callScore['sigma']*(cluster$x2-cluster$x1)/5e6
    return(callScore)
  })
  return(t(ret))
}

callConsistencies = function(clusters, shift=0, N=1000) {
  ret = lapply(1:nrow(clusters), function(row) {
    superFreq:::clusterToConsistencies(clusters[row,], shift=shift, N=N)
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


findShift = function(clusters, nShifts=100, maxShift = 3, cpus=1, plot=T, plotPloidy=F, ploidyPrior=NULL) {
  #set up shifts more densely around 0
  nShiftsHalf = round(nShifts/2)
  shifts = ((0:nShiftsHalf)/nShiftsHalf)^2*maxShift
  shifts = c(-rev(shifts), shifts[-1])

  #get the scores from consistency for each shift
  shiftConsistencies = mclapply(shifts, function(shift) {
    superFreq:::callConsistencies(clusters, shift=shift, N=1000)
  }, mc.cores=cpus)

  #find the shift with best score
  shiftScores = sapply(shiftConsistencies, superFreq:::scoreConsistencies)
  if ( !is.null(ploidyPrior) ) {
    basePloidy = 2*sum((clusters$x2-clusters$x1)*2^clusters$M)/sum((clusters$x2-clusters$x1))
    ploidies = basePloidy*2^shifts
    #add penalty to the discrepancy proportional to square of distance from prior.
    #but don't go too far above the max discrepancy before, for viz reasons
    shiftScores = pmin(max(shiftScores)+0.1, shiftScores + (ploidies-ploidyPrior)^2)
  }
  bestShift = shifts[shiftScores == min(shiftScores)][1]

  #return if no shift was best
  if ( bestShift == 0 ) {
    if ( plot ) superFreq:::plotShiftedMAFLFC(clusters)
    if ( plotPloidy ) {
      basePloidy = 2*sum((clusters$x2-clusters$x1)*2^clusters$M)/sum((clusters$x2-clusters$x1))
      ploidies = basePloidy*2^shifts
      main = 'ploidy fit'
      if ( !is.null(ploidyPrior) ) main = paste0('ploidy fit (with prior=', ploidyPrior, ')')
      superFreq:::plotColourScatter(ploidies, shiftScores, cex=2, xlim=basePloidy*c(0.2, 3), xlab='ploidy', ylab='fit discrepancy', main='ploidy fit')
      lines(ploidies, shiftScores, lwd=2, col=mcri('blue'))
      ymax = max(shiftScores)
      yscale = ymax - min(shiftScores)
      segments(basePloidy, 0, basePloidy, ymax-yscale*0.05, lwd=1.5, col=mcri('orange'))
      text(basePloidy, ymax, round(basePloidy, digits=2), font=2, col=mcri('orange'))
    }
    return(clusters)
  }

  #if a non-zero shift is better scored, renormalise and rerun to better find local minimum
  #with the denser shifts around 0.
  clusters$M = clusters$M + bestShift
  clusters = findShift(clusters, nShifts=nShifts, maxShift=maxShift, cpus=cpus, plot=plot, plotPloidy=plotPloidy, ploidyPrior=ploidyPrior)
  return(clusters)
}

plotShiftedMAFLFC = function(clusters, shift=0) {
  clusters$M = clusters$M + shift
  return(plotMAFLFC(clusters))
}



findShiftManually = function(clusters, nShifts=100, maxShift = 3, cpus=1, plot=T, plotPloidy=F) {
  #set up shifts more densely around 0
  nShiftsHalf = round(nShifts/2)
  nullShifts = ((0:nShiftsHalf)/nShiftsHalf)^2*maxShift
  nullShifts = c(-rev(nullShifts), nullShifts[-1])

  #ask user for the shift
  userHappy = F
  userShift = 0
  while( !userHappy ) {
    #calculate state for user shift
    cat('Calculating inconsistency scores')
    shifts = nullShifts + userShift
    shiftConsistencies = mclapply(shifts, function(shift) {
      cat('.')
      callConsistencies(clusters, shift=shift, N=1000)
    }, mc.cores=cpus)
    shiftScores = sapply(shiftConsistencies, scoreConsistencies)
    cat('done.\n')
    
    #plot info about current shift
    resetMargins()
    layout(matrix(c(1, 1, 2, 3), nrow=1))
    plotShiftedMAFLFC(clusters, shift=userShift)
    plotCR(clusters)
    plot(shifts, shiftScores, type='c', xlab='LFC shift', ylab='inconsistency score')
    text(shifts, shiftScores, 1:length(shifts))
    ploidy = sum(2*2^(clusters$M + userShift)*(clusters$x2-clusters$x1))/sum(clusters$x2-clusters$x1)
    cat('Ploidy with shift ', userShift, ' is ', round(ploidy, 3), '.\n', sep='')
    cat('Inconsistency score is ', shiftScores[21], '.\n', sep='')
    cat('Minimum score is ', min(shiftScores), ' at index ', which(shiftScores==min(shiftScores)),
        ', with shift ', round(shifts[which(shiftScores==min(shiftScores))], 3), '.\n', sep='')

    #ask for new shift, return if user happy.
    i = -1
    while ( !(i %in% 0:length(shifts)) ) {
      i = as.numeric(readline('type index of new shift from the right plot. 0 --> use this normalisation.\nindex: '))
    }
    if ( i == 0 ) {
      i = 21
      userHappy = T
    }
    userShift = shifts[i]
  }
  cat('Final ploidy with shift ', userShift, ' is ', round(ploidy, 3), '.\n', sep='')
  layout(1)
  
  bestShift = userShift

  #return if no shift was best
  if ( bestShift == 0 ) {
    if ( plot ) plotShiftedMAFLFC(clusters)
    return(clusters)
  }

  #if a non-zero shift is better scored, renormalise and rerun to better find local minimum
  #with the denser shifts around 0.
  clusters$M = clusters$M + bestShift
  return(clusters)
}
