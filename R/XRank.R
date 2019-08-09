#########################################################################
#
# The main function GLFC(fit, pval=0.1, mode='t', priors='empirical'
#                        coefs=0, plot=T)
# takes the output from eBayes, and outputs a matrix of guaranteed log
# fold changes, at the confidence of the p-value, for the coefficients
# specified. The log fold change can have (mode) 't' or 'normal' distributed
# likelyhoods, and the priors can be determined 'empirical', 'flat', or
# supplied diectly, for example from output from 'getPriors'.
# Plot determines whether a volcano plot for the first
# given coefficient is shown or not.
#
# plotVolcano(glfc, fit, coef, main = 'volcano plot', print = 0)
# takes a list of guaranteed log fold changes (so a column of the output
# from GLFC), and the corresponding fit and coefficient, and shows
# how the log fold changes are modified.
# 'print' says how many of the top (and bottom) glfc are to be labeled
# in the plot with the gene names in the fit object.
#
# Last changed 3 July 2013
#
#########################################################################

#The main function, to be called with a fit object
#created through the limma pipeline
GLFC = function(fit, pval = 0.1, mode='t', priors='empirical', FDR = F,
                plot=T, quiet = F, bestGuess='default') {
  coefs = colnames(fit$coefficients)


  if ( pval ==1 && bestGuess == 'default' ) {
    bestGuess = T
    if ( !quiet ) cat('Defaulting bestGuess to True, as pval is 1.\n')
  }
  else if ( bestGuess == 'default' ) bestGuess = F
  else if ( pval != 1 && bestGuess == T )
    cat('\nWARNING: bestGuess only recommended for pval = 1. Interpret with caution!\n\n')
  
  #get the priors
  if ( is.character(priors) && priors == 'empirical' ) {
    #call the emprical prior function.
    if ( !quiet ) cat('Preparing empirical priors...')
    priors = getPriors(fit, coefs = coefs, plot=plot)
    if ( !quiet ) cat('done.\n')
  }
  else if ( is.character(priors) && priors == 'flat' ) {
    #call the flat prior function.
    max = sapply(1:length(coefs), function(i) 1.1*max(abs(fit$coefficients[,i])))
    priors = lapply( 1:length(coefs), function(i) getFlatPrior(max = max[i], plot=plot))
  }
  else if ( is.list(priors) ) {
    #if user supplied prior, do a sanity check.
    if ( length(priors) != length(coefs) ) {
      cat(length(priors),' priors, and ', length(coefs),
          ' coefficients. Need matching numbers.\n',sep='')
      return(-1)
    }
    else if ( !quiet ) cat('Using provided priors.\n')
  }
  else {
    cat('Could not make sense of prior input. Please use \'empirical\', \'flat\' or provide a list of prepared priors.\n')
    return(-1)
  }
  
  if ( !quiet ) cat('Calculating GLFC...')

  #calculate binsizes once for all
  getBinsize = function(prior) {
    nbins = length(prior$mids)
    prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }
  binsize = lapply(priors, getBinsize)
  names(binsize) = names(priors)

  #loop through the genes and coefficients, calculating the safe
  #ratio for each by calling the GR function.
  #The uncertainty in LFC is taken from the t-statistics or the
  #s2.post and stdev.unscaled of the fit object, depending on the input option.
  getGLFC = function(coef, gene) {
    GR(2^fit$coefficients[gene,coef], p=pval, prior = priors[[coef]],
           mode=mode, df = fit$df.total[gene], t=abs(fit$t[gene,coef]),
           d = sqrt(fit$s2.post[gene])*fit$stdev.unscaled[gene, coef],
           bestGuess=bestGuess, binsize = binsize[[coef]])
  }
  
  getColumn = function(coef) {
    if ( !quiet ) cat(coef, '..', sep='')
      do.call(rbind,
              lapply(1:length(rownames(fit)), function(gene)
                     getGLFC(coef, gene)))
  }

  data = do.call(cbind,
    lapply(coefs,
           getColumn))

  GLFCs = data[,(1:length(coefs))*2-1]
  ps = data[,(1:length(coefs))*2]
  
  if ( !quiet ) cat('done.\n')

  #name the columns, for prettier output.
  if ( length(coefs) > 1 ) colnames(GLFCs) = coefs

  #plot the volcano plot
  
  if ( plot ) {
    if ( length(coefs) > 1 )
      plotVolcano(GLFCs[,2], fit, coefs[2], print = 5,
                  main = coefs[2])
    else
      plotVolcano(GLFCs, fit, coefs, print = 5,
                  main = colnames(fit$coefficients)[coefs[1]])
  }

  #Change from FC to LFC
  GLFCs = log2(GLFCs)
  #name the rows.
  if (  length(coefs) == 1  )
    names(GLFCs) = rownames(fit)
  else
    rownames(GLFCs) = rownames(fit)

  return(list(glfc=GLFCs, p=ps))
}

#This function calculates the prior from the data supplied in the fit object.
#The bins are denser around LFC = 0, so default number of bins should be enough.
getPriors = function(fit, coefs = 0, nbins = 50, plot=T, mode='empirical') {
  #if no coefficients specified, get priors for all.
  if ( length(coefs) == 1 && coefs == 0 ) coefs = colnames(fit$coefficients)
  N = length(coefs)

  #find maximum LFC, to know the range to take data over.
  range = lapply(coefs, function(coef) 1.1*max(abs(fit$coefficients[,coef])))
  names(range) = coefs

  #Set the breaks for the bins in the histogram.
  #Smaller bins close to 0, where accuracy is more important,
  #and density is higher
  breaks = lapply(coefs, function(coef)
    (((-nbins):(nbins))/nbins)^2*sign((-nbins):(nbins))*range[[coef]] )
  names(breaks) = coefs

  #Now get the histograms, which will be the priors
  priors = lapply(coefs, function(coef) {
    if ( mode == 'posterior' & 'best.guess' %in% names(fit) ) hist(fit$best.guess[,coef], breaks = breaks[[coef]], plot=plot)
    else hist(fit$coefficients[,coef], breaks = breaks[[coef]], plot=plot)
  })
  names(priors) = coefs

  #Return the priors
  return(priors)
}

#A function generating some data evenly distributed in LFC up to
#the supplied maximum LFC, and then transforms it into a flat prior.
getFlatPrior = function(max = 15, nbins = 100, plot=T) {
  #Set the breaks for the bins in the histogram.
  #Smaller bins close to 0, where accuracy is more important,
  #and density is higher
  breaks = (((-nbins):(nbins))/nbins)^2*
    sign((-nbins):(nbins))*max

  #fake a flat dataset
  flat = runif(100000, -max, max)

  #Now get the histogram, which will be the prior
  prior = hist(flat, breaks = breaks, plot=plot)

  #Return the priors
  return(prior)
}

#This is the heart of the algorithm, calculating the guaranteed ratio
#for a given p-value, measured ratio, a given prior,
#and a given uncertainty in the LFC. The last may be normally or t distributed.
#Further there are options to plot the prior, prior probability distribution
#and posterior probability distribution, with shaded area representing the
#requested p-value. This plot is very helpful in understanding the algorithm
#when a single gene is studied, but cumbersome when large sets of genes are analysed.
GR = function(r, d=0.1, p=0.1, prior = 'flat', mode='normal',
  df=4, t=1, verbose=F, zoom = F, bestGuess = F, binsize = 0) {
  #just a shorthand notation.
  v = verbose

  #get the flat prior if needed.
  if ( is.character(prior) &&  prior == 'flat' ) {
    max = 2*(abs(log2(r))+d+0.5)
    prior = getFlatPrior(max = max, plot=F)
  }
  
  #set the grid for the prior.
  x = prior$mids
  #The density, not depending on binwidth.
  priorD = prior$density

  #calcualte binsize on the grid
  if ( !is.list(binsize) ) {
    nbins = length(x)
    binsize = prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }

  #set the prior probability distribution of the LFC on the grid.
  #Both the density (useful for multiplying probabilities) and
  #"counts" (useful for integrating) are kept track of.
  if ( mode == 'normal' ) {
    flatD = dnorm(x, log2(r), d)
  }
  else if ( mode == 't' ) {
    flatD = dt((x-log2(r))*t/abs(log2(r)), df, 0)
  }

  #find the normalised (over all bins) posterior prob dist
  #Multiply the densities to get the posterior probability density.
  postD = flatD*priorD
  #Transformt o counts for easier integration.
  postC = postD*binsize
  #Normalise to give probabilities.
  postCNorm = postC/sum(postC)

  #Find the LFC where the integrated probability of the posterior is p/2.
  if ( r > 1 ) x0 = findGLFCandPV(postCNorm, prior$breaks, p/2)
  else x0 = findGLFCandPV(postCNorm, prior$breaks, 1 - p/2)

  #split apart the 'p-value'
  p = x0[2]
  x0 = x0[1]

  if (v) plotShrunk(r, d=d, shade=2^x0, zoom = zoom, prior = prior,
                    mode=mode, df=df, t=t, side='left', reversed = F)

  if ( sign(x0*log(r)) == -1 && !bestGuess ) return(c(1, p))

#  if ( reversed ) x0 = -x0
  return(c(2^x0, p))
}

#This function plots the p-values (x) against the LFCs (y).
#The largest 'print' positive and negative GLFCs can be marked by their gene
#names in the fit object if desired.
plotVolcano = function(fit, coef=1, print = 20, main = NA, shrink=T, line='nGenes',
  xlab='corrected LFC', ylab='-log10(p-value)', specialGenes=c(), col='red', ...) {
  #get the name of the column if only index provided
  if ( is.numeric(coef) ) coef = colnames(fit)[coef]
  print=min(print, nrow(fit))
  #get the FCs and p values
  bgs = fit$best.guess[,coef]
  cos = fit$coefficients[,coef]
  ps = fit$p.value[, coef]

  #decide for a range on the x-axis
  xmax = max(1, 1.5*max(abs(bgs)))
  xlim = c(-xmax, xmax)

  if ( is.na(main) ) main = coef

  ps = pmax(ps, pmin(min(ps[ps > 0], 1e-10)))
  
  #plot the points and segments
  plot(1, type='n', xlim = xlim, ylim=c(0, max(-log10(ps))), xlab=xlab, ylab=ylab, main = main, ...)
  if ( line == 'nGenes' ) segments( 2*xlim[1], log10(nrow(fit)), 2*xlim[2], log10(nrow(fit)), cex=1, col='orange', lty=2)
  if ( line == 'fdr' ) {
    fdr = p.adjust(ps, method='fdr')
    cut = -log10(max(ps[fdr <= 0.05]))
    segments(2*xlim[1], cut, 2*xlim[2], cut, cex=1, col='orange', lty=2)
  }
  segments(bgs, -log10(ps), cos, -log10(ps), cex=0.4, col=rgb(0,0,0,0.15))
  points(cos, -log10(ps), pch=16, cex=0.5)
  points(bgs, -log10(ps), pch=16, cex=0.7, col=col)

  legend('bottomright', c('measured', 'best guess'), pch=16, pt.cex=c(0.5, 0.8), col=c('black', 'red'))

  #add labels for top DE tags.
  if ( print > 0 ) {
    names = rownames(fit)
    tops = which(cos > 0)[order(fit$XRank[cos > 0,coef])[1:print]]
    bots = which(cos < 0)[order(fit$XRank[cos < 0,coef])[1:print]]

    ymax = -log10(min(ps))
    if ( shrink ) cex = 0.6 + 0.4*(print:1)/print
    else cex = 0.8
    shift = (0.5+1.5*cex)*ymax/100
    text(bgs[tops], -log10(ps[tops])+shift, names[tops], col='red', cex=cex)
    text(bgs[bots], -log10(ps[bots])+shift, names[bots], col='red', cex=cex)
    special = which(rownames(fit) %in% specialGenes)
    if ( length(special) > 0 ) {
      points(bgs[special], -log10(ps[special]), col=rgb(0, 0.8, 1), cex=1, pch=16)
      text(bgs[special], -log10(ps[special])+ymax/50, names[special], col=rgb(0, 0.8, 1), cex=1)
    }
  }
}


#The plot function for the probability distributions.
#Usually better to call GR(..., plot = TRUE) than calling this directly.
plotShrunk = function(r, d=0.1, prior = 'empirical', coef=1, shade=0, mode='normal',
  df=4, t=1, side='left', zoom = F, reversed = F) {

  if ( reversed ) {
    r = 1/r
    shade = 1/shade
    if ( side == 'left' ) side = 'right'
    else side = 'left'    
  }
  
  #get the priors
  if ( is.character(prior) && prior == 'empirical' ) {
    cat('Preparing empirical priors... ')
    prior = getPriors(fit, coefs = coef, plot=F)
    cat('done.\n')
  }
  else if ( is.character(prior) && prior == 'flat' ) {
    max = 2*(abs(log2(r))+d+0.5)
    prior = getFlatPrior(max = max, plot=F)
  }
  else cat('Using provided prior.\n')


  #set the grid on lFC
  x = prior$mids

  #define the prior distribution and counts per bin.
  priorD = prior$density

  #normalise the density to 1.
  priorDNorm=priorD/sum(priorD)

  #calculate P(true lFC | observed lFC) on a flat prior
  #using the provided observed lFC and variance, and normal or t dist.
  if ( mode == 'normal' ) {
    flatD = dnorm(x, log2(r), d)
    w = d
  }
  else if ( mode == 't' ) {
    flatD = dt((x-log2(r))*t/abs(log2(r)), df, 0)
    w = abs(log2(r))/t
  }
  else {
    cat('Set mode to \'normal\' or \'t\'\n.')
    return(-1)
  }

  #Normalise and find the posterior P(true lFC | observed lFC)
  #on the given prior.
  flatDNorm = flatD/sum(flatD)
  postD = flatD*priorD
  postDNorm = postD/sum(postD)

  #zoom in on the lower part of the plot if so desired
  if ( zoom ) {
    limit = max(priorDNorm, flatDNorm, postDNorm)/20
    zoomfnc = function(x) ifelse(x<limit, x*10, x+limit*9)
    priorDNorm = zoomfnc(priorDNorm)
    flatDNorm = zoomfnc(flatDNorm)
    postDNorm = zoomfnc(postDNorm)
  }

  #plot the prior, flat P(true|obs) and posterior P(true, obs)
  plot(2^x, priorDNorm, type='l', log='x', xlim=c(min(0.8, (r*2^(-w))/1.4), max(1.2,(r*2^(w))*1.4)),
       ylim=c(0.001, max(priorDNorm, flatDNorm, postDNorm)))
  lines(2^(x), flatDNorm, col='blue')
  lines(2^(x), postDNorm, col='red')

  #shade the zoomed in area
  if ( zoom ) polygon(c(0.00001, 2^(x), 100), c(0, x/x*limit*10, 0), col = rgb(0, 0.8, 0, 0.1))

  #shade part of the area under the posterior.
  if ( shade != 0 & side == 'left')
    polygon(2^(x), postDNorm*as.integer(2^(x) < shade), col = rgb(0.7, 0, 0, 0.5))
  if ( shade != 0 & side == 'right')
    polygon(2^(x), postDNorm*as.integer(2^(x) > shade), col = rgb(0.7, 0, 0, 0.5))
}


findGLFC = function(dist, x, p) {
  sums = cumsum(dist)
  
  i = max(which(sums < p))

  x0 = x[i] + (p - sums[i])*(x[i+1] - x[i])/dist[i+1]

  return(x0)
}

findGLFCandPV = function(dist, breaks, p) {
  dist = as.numeric(dist)
  breaks = as.numeric(breaks)
  #integrated probability distribution
  sums = cumsum(as.numeric(dist))

  #the break that bring the integral over p
  i = max(which(sums < p))

  #interpolate to find the x between the bins best fitting the limit.
  x0 = breaks[i+1] + (p - sums[i])*(breaks[i+2] - breaks[i+1])/dist[i+1]

  #integrated probability to be below 0
  p0 = interpolate(c(-Inf, breaks), sums, 0)
  #find smallest side (above or below 0).
  #multiply by 2 for double side test.
  p0 = 2*min(p0, 1-p0)

  return(c(x0, p0))
}








###########################################################################
#
# Posterior analysis
#
##########################################################################


getPosterior = function(r, d=0.1, prior = 'flat', mode='normal',
  df=4, t=1, verbose=F, zoom = F, binsize = 0) {
  #just a shorthand notation.
  v = verbose

  #get the flat prior if needed.
  if ( is.character(prior) &&  prior == 'flat' ) {
    max = 2*(abs(log2(r))+d+0.5)
    prior = getFlatPrior(max = max, plot=F)
  }
  
  #set the grid for the prior.
  x = prior$mids
  #The density, not depending on binwidth.
  priorD = prior$density

  #calcualte binsize on the grid
  if ( !is.list(binsize) ) {
    nbins = length(x)
    binsize = prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }

  #set the prior probability distribution of the LFC on the grid.
  #Both the density (useful for multiplying probabilities) and
  #"counts" (useful for integrating) are kept track of.
  if ( mode == 'normal' ) {
    flatD = dnorm(x, log2(r), d)
  }
  else if ( mode == 't' ) {
    flatD = dt((x-log2(r))*t/abs(log2(r)), df, 0)
  }

  #find the normalised (over all bins) posterior prob dist
  #Multiply the densities to get the posterior probability density.
  postD = flatD*priorD
  #Transformt to counts for easier integration.
  postC = postD*binsize
  #Normalise to give probabilities.
  postCNorm = postC/sum(postC)

  return(postCNorm)
}

posteriors = function(fit, mode='t', priors='empirical', FDR = F,
                     coefs = 0, plot=T, quiet = F, cpus=1) {
  if ( coefs[1] == 0 ) coefs = 1:ncol(fit)
  if ( is.numeric(coefs) ) coefs = colnames(fit)[coefs]

  #get the priors
  if ( is.character(priors) && priors == 'empirical' ) {
    #call the emprical prior function.
    if ( !quiet ) cat('Preparing empirical priors... ')
    priors = getPriors(fit, coefs = coefs, plot=plot, mode='empirical')
    if ( !quiet ) cat('done.\n')
  }
  else if ( is.character(priors) && priors == 'posterior' ) {
    #call the emprical prior function.
    if ( !quiet ) cat('Preparing empirical priors... ')
    priors = getPriors(fit, coefs = coefs, plot=plot, mode='posterior')
    if ( !quiet ) cat('done.\n')
  }
  else if ( is.character(priors) && priors == 'flat' ) {
    #call the flat prior function.
    max = sapply(1:length(coefs), function(i) 1.1*max(abs(fit$coefficients[,i])))
    priors = lapply( 1:length(coefs), function(i) getFlatPrior(max = max[i], plot=plot))
  }
  else if ( is.list(priors) ) {
    #if user supplied prior, do a sanity check.
    if ( length(priors) != length(coefs) ) {
      cat(length(priors),' priors, and ', length(coefs),
          ' coefficients. Need matching numbers.\n',sep='')
      return(-1)
    }
    else if ( !quiet ) cat('Using provided priors.\n')
  }
  else {
    cat('Could not make sense of prior input. Please use \'empirical\', \'flat\' or provide a list of prepared priors.\n')
    return(-1)
  }
  
  if ( !quiet ) cat('Calculating posteriors:\n')

  #calculate binsizes once for all
  getBinsize = function(prior) {
    nbins = length(prior$mids)
    prior$breaks[2:(nbins+1)] - prior$breaks[1:nbins]
  }
  binsize = lapply(priors, getBinsize) 

  #loop through the genes and coefficients, calculating the safe
  #ratio for each by calling the GR function.
  #The uncertainty in LFC is taken from the t-statistics or the
  #s2.post and stdev.unscaled of the fit object, depending on the input option.
  getPosteriors = function(rep, gene) {
    getPosterior(2^fit$coefficients[gene, rep], prior = priors[[rep]],
                 mode=mode, df = fit$df.total[gene], t=abs(fit$t[gene,rep]),
                 d = sqrt(fit$s2.post[gene])*fit$stdev.unscaled[gene, rep],
                 binsize = binsize[[rep]])
  }
  
  getMatrix = function(rep) {
    if ( !quiet ) cat(rep, '...', sep='')
    do.call(rbind,
            mclapply(1:length(rownames(fit)), function(gene) getPosteriors(rep, gene), mc.cores=cpus))
  }

  data = lapply(coefs, getMatrix)
  
  if ( !quiet ) cat('done.\n')

  #name the columns, for prettier output.
  names(data) = coefs
  for ( coef in coefs ) {
    rownames(data[[coef]]) = rownames(fit)
    colnames(data[[coef]]) = priors[[coef]]$mids
  }

  return(list(data = data, prior=priors, fit = fit))
}

generateNullFit = function(fit) {
  ret = fit
  ret$coefficients = fit$coefficients/fit$t*rt(length(fit$coefficients), fit$df.total, 0)
  return(ret)
}

interpolates = function(x, y, x0s) {
  y0s = sapply(x0s, function(x0) interpolate(x, y, x0))
  return(y0s)
}

interpolate = function(x, y, x0) {
  if ( x0 < min(x) ) return(y[1])
  if ( x0 > max(x) ) return(y[length(y)])
  if (x0 %in% x) return(y[which(x == x0)])
  
  lowI = max(which(x < x0))
  highI = min(which(x > x0))
  lowX = x[lowI]
  highX = x[highI]
  lowY = y[lowI]
  highY = y[highI]
  ret = lowY + (x0-lowX)/(highX-lowX)*(highY-lowY)
  return(ret)
}

#finds and returns the difference from the prior distribution expected, given the
#t-distribution, variances in fit, and all null hypothesis true.
getPriorNullDeviation = function(fit, xbreaks, wbreaks=10) {

  return(
    lapply(1:ncol(fit$coefficients), function(coef)
           getDeviation(fit, xbreaks, coef, wbreaks=wbreaks)
           )
    )
}

getDeviation = function(fit, xbreaks, coef, wbreaks=10) {
  nbins = length(xbreaks)-1
  xs = (xbreaks[1:nbins] + xbreaks[2:(nbins+1)])/2
  binsize = xbreaks[2:(nbins+1)] - xbreaks[1:nbins]
  
  ws = fit$coefficients[,coef]/fit$t[,coef]
  #breaks = unique(quantile(ws, seq(0, 1, 1/wbreaks)))
  breaks = getWBreaks(ws, wbreaks)
  
  wHist = hist(ws, breaks=breaks)
  wBin = sapply(ws, function(w)
    max(which(wHist$breaks <= w))
    )
  wList = lapply(1:length(wHist$mids), function(i)
    which(wBin == i)
    )
  
  expectedMx =
    do.call(rbind,
            lapply(wList, function(is)
                   binsize*sapply(xs, function(x)
                                  sum(dt(x/ws[is], fit$df.total[1], 0)/ws[is])
                                  )
                   )
            )
  
  rownames(expectedMx) = wHist$mids
  colnames(expectedMx) = xs
  
  return(expectedMx)
}

getWBreaks = function(ws, breaks=10) {
  return(unique(quantile(ws, sort(c(seq(0, 1, 1/breaks), 1 - 1/breaks/2, 1 - 1/breaks/4)))))
}


DEposterior = function(fit, coef, wbreaks=10, FDRcut = 0.5) {
  priors = getPriors(fit, plot=F)
  xbreaks = priors[[coef]]$breaks
  minC = min(fit$coefficients[,coef])
  maxC = max(fit$coefficients[,coef])
  nx = 100
  x = (xbreaks[1:(length(xbreaks)-1)] + xbreaks[2:length(xbreaks)])/2
  null = getPriorNullDeviation(fit, xbreaks, wbreaks=wbreaks)[[coef]]

  ws = fit$coefficients[,coef]/fit$t[,coef]
  wbreaks = getWBreaks(ws, wbreaks)
  
  wHist = hist(ws, breaks=wbreaks)
  wBin = sapply(ws, function(w)
    max(which(wHist$breaks <= w))
    )
  wList = lapply(1:length(wHist$mids), function(i)
    which(wBin == i)
    )

  DEpost = nullRatio = meas = correction = matrix(0, ncol=ncol(null), nrow=nrow(null))
  for ( i in 1:length(wList) ) {
    if ( length(wList[[i]]) == 0 ) next
    meas[i,] = hist(fit$coefficients[wList[[i]],coef], plot=F,
      breaks=xbreaks)$counts

    #ratio of true null hypothesis
    Nn = scalarNorm(multiBlur(meas[i,], range=2, repeats=5),
      multiBlur(null[i,], range=2, repeats=5))
        
    DEpost[i,] = noneg(multiBlur(meas[i,], range=1, repeats=3) - Nn*(multiBlur(null[i,]+0*sqrt(null[i,]), range=1, repeats=3)))

    TDR = function(i, j) sum(DEpost[i,abs(x) >= abs(x[j])])/
      sum(DEpost[i,abs(x) >= abs(x[j])] + null[i,abs(x) >= abs(x[j])])

    correction[i,] = sapply(1:length(x), function(j) TDR(i, j))
    DEpost[i,] = DEpost[i,]*ifelse(correction[i,] < 1 - FDRcut, 0, 1)

    cat('Ratio of DE genes at width around ',signif(wHist$mids[i], digits=2), ' is ',
        round(100*sum(DEpost[i,])/sum(meas[i,])), '%\n', sep='')
    
    nullRatio[i,] = Nn*null[i,]/meas[i,]
  }

  cat('total number of DE genes is estimated to ',round(sum(DEpost)),
      ', which is ', round(100*sum(DEpost)/length(wBin)), '% of all genes.\n', sep='')

  colnames(DEpost) = colnames(null) = colnames(correction) = colnames(meas) = x
  rownames(DEpost) = rownames(null) = rownames(correction) = rownames(meas) = wHist$mids

  plot(x, colSums(null), type='l', xlim=c(-1.5,1.5),
       col=mcri('blue'), lwd=2)
  lines(x, colSums(meas), col=mcri('red'), lwd=2)
  lines(x, colSums(DEpost), col=mcri('green'), lwd=2)
  legend('topleft', c('null hypothesis', 'measured', 'DE'),
         col=c(mcri('blue'), mcri('red'), mcri('green')), lwd=3)
  
  return(list(DE = DEpost, null = null, measured = meas, TDR = correction))
}

adjustBG = function(BG, fit, dep, coef=2) {
  wBin = function(w) max(c(1, which(as.numeric(rownames(dep$DE)) <= w)))
  xBin = function(x) max(c(1, which(as.numeric(colnames(dep$DE)) <= x)))

  correction = sapply(1:length(BG), function(i) {
    w = wBin(fit$coefficients[i,coef]/fit$t[i,coef])
    x = xBin(fit$coefficients[i,coef])
    return(1/(dep$null[w, x]/dep$DE[w,x] + 1))
    })

  return(BG*correction)
}


postRank = function(posts, quiet=F, cpus=1) {
  if ( !quiet ) cat('Calculating expected ranks...')

  #set up return matrix
  retRank = matrix(0, ncol=length(posts$data), nrow=nrow(posts$data[[1]]))
  retMeanRank = matrix(0, ncol=length(posts$data), nrow=nrow(posts$data[[1]]))
  colnames(retRank) = colnames(retMeanRank) = names(posts$data)
  rownames(retRank) = rownames(retMeanRank) = rownames(posts$data[[1]])

  for ( coef in names(posts$data) ) {
    if ( !quiet ) cat(coef, '..', sep='')
    mx = posts$data[[coef]]
    mx = mx[,1:50] + mx[,100:51]
    xs = posts$prior[[coef]]$mids
    xs = xs[1:50]
    F = colSums(mx)

    ranks = cumsum(F) - F/2 #TODO check that this added -F/2 works in practice!!
    binsize = F
    cumPost = do.call(rbind, mclapply(1:nrow(mx), function(gene) cumsum(mx[gene,]), mc.cores=cpus))
    
    currentRank = 1
    ranking = rep(1, nrow(mx))
    for (i in 1:length(ranks)) {
      if ( floor(ranks[i]) < currentRank ) next
      n = floor(ranks[i]) - currentRank + 1
      tops = order(-cumPost[,i])[1:n]
      ranking[currentRank:(currentRank+n-1)] = tops
      cumPost[tops,] = cumPost[tops,]*0
      currentRank = currentRank + n
    }
    
    meanRank = sapply(1:nrow(mx), function(gene) sum(ranks*mx[gene,]))
    if ( T | coef == names(posts$data)[1] ) {
      retRank[,coef] = ranking
      retMeanRank[,coef] = meanRank
      #retRank = ranking
      #retMeanRank = meanRank
    }
    else {
      retRank = cbind(retRank, ranking)
      retMeanRank = cbind(retMeanRank, meanRank)
    }
  }

#  if ( length(posts$data) > 1 ) {
#    colnames(retRank) = colnames(retMeanRank) = names(posts$data)
#    rownames(retRank) = rownames(retMeanRank) = rownames(posts$data[[1]])
#  }
#  else {
#    names(retRank) = names(retMeanRank) = rownames(posts$data[[1]])
#  }

  if ( !quiet ) cat('done.\n')
  return(list(ranking = retRank, XRank = retMeanRank))
}

plotPost = function(posts, i, truth=F, BG=F) {
  p = posts$data[[1]][,i]
  xs = posts$prior[[1]]$mids

  plot(xs, multiBlur(p, 1, 2), type='l', xlim=c(-1.5, 1.5))
  if ( length(truth) > 1 ) {
    segments( truth[i], 0, truth[i], 1, col=mcri('red'), lwd=2)
    cat('true ranking is ', rank(-abs(truth))[i], '.\n', sep='')
  }
  segments(fit$coefficients[i, 2], 0, fit$coefficients[i, 2], 1, mcri('blue'), lwd=2)
  if ( length(BG) > 1 )
  segments(BG$glfc[i], 0, BG$glfc[i], 1, col=mcri('green'), lwd=2)

  legend('topleft', c('meas', 'truth', 'BG'), col=c(mcri('blue'), mcri('red'), mcri('green')), lwd=2)
}

#' Turns colours into similar colours from the Murdoch Childrens Research Institute palette.
#'
#' @details The function return a colour from the MCRI palette if match found, otherwise the input is returned unchanged. call mcri() to see the available colours.
#'
#' @param col character or numeric. Most common colours such as "green" or 'blue' are converted to MCRI. Numbers return the colour matching that number.
#' @param al The alpha parameter: the opaqueness. Numeric between 0 and 1.
#'
#' @export
#' @examples
#' plot(1:9, rep(1,9), pch=16, cex=10, col=mcri(0:8))
#' text(1:9, rep(1.1,9), c('blue', 'orange', 'green', 'magenta',
#'                         'cyan', 'red', 'violet', 'darkblue', 'darkred'),
#'      col=mcri(0:8), cex=1.5, font=2)
#' 
mcri = function(col='deafult', al=1) {
  if ( length(col) == 0 | length(al) == 0 ) return(character(0))
  if ( length(col) == 1 && col[1] == 'deafult' ) {
    cat('Use: mcri(\'colour\'), returning an official MCRI colour.\nAvailable MCRI colours are:\n\ndarkblue\nblue\nlightblue\nazure\ngreen\norange\nviolet\ncyan\nred\ndarkred\nmagenta (aka rose).\n\nReturning default blue.\n')
    return(mcri('blue'))
  }
  if ( length(col) > 1 & length(al) == 1 ) return(sapply(col, function(c) mcri(c, al)))
  else if ( length(col) > 1 & length(al) == length(col) ) return(sapply(1:length(col), function(i) mcri(col[i], al[i])))
  else if ( length(col) > 1 ) {
    warning('length of col and al mismatch in mcri()')
    return(sapply(col, function(c) mcri(c, al[1])))
  }

  if ( length(col) == 1 & length(al) > 1 ) return(sapply(al, function(alpha) mcri(col, alpha)))
  
  if ( is.numeric(col) ) {
    col = (col %% 9) + 1
    if ( col == 1 ) col = 'blue'
    else if ( col == 2 ) col = 'red'
    else if ( col == 3 ) col = 'green'
    else if ( col == 4 ) col = 'magenta'
    else if ( col == 5 ) col = 'cyan'
    else if ( col == 6 ) col = 'orange'
    else if ( col == 7 ) col = 'violet'
    else if ( col == 8 ) col = 'darkblue'
    else if ( col == 9 ) col = 'darkred'
    else col = 'black'
  }
  ret = 0
  if ( col == 'darkblue') ret = rgb(9/255, 47/255, 94/255, al)
  if ( col == 'blue') ret = rgb(0, 83/255, 161/255, al)
  if ( col == 'lightblue') ret = rgb(0, 165/255, 210/255, al)
  if ( col == 'azure') ret = rgb(0, 173/255, 239/255, al)
  if ( col == 'green') ret = rgb(141/255, 198/255, 63/255, al)
  if ( col == 'orange') ret = rgb(244/255, 121/255, 32/255, al)  
  if ( col == 'violet') ret = rgb(122/255, 82/255, 199/255, al)  
  if ( col == 'cyan') ret = rgb(0/255, 183/255, 198/255, al)  
  if ( col == 'red') ret = rgb(192/255, 80/255, 77/255, al)  
  if ( col == 'darkred') ret = rgb(96/255, 40/255, 38/255, al)  
  if ( col == 'magenta' | col == 'rose') ret = rgb(236/255, 0/255, 140/255, al)
  if (ret == 0 ) ret = do.call(rgb, as.list(c(col2rgb(col)/255, al)))
  return(ret)
}


bestGuess = function(fit, coefs = 0, quiet=F, cpus=1) {
  if ( !quiet ) cat('Calculating best guess...')
  fit$best.guess = do.call(cbind,
    lapply(coefs, function(coef) {
      if ( !quiet ) cat(coef, '..', sep='')
      unlist(mclapply(1:nrow(fit), function(gene)
                      findGLFCandPV(fit$posterior[[coef]][gene,],
                                    fit$prior[[coef]]$breaks, 0.5)[1], mc.cores=cpus))
    }
           )
    )
  colnames(fit$best.guess) = names(fit$posterior)
  rownames(fit$best.guess) = rownames(fit)
  if ( !quiet ) cat('done.\n')
  return(fit)
}

XRank = function(fit, coefs = 0, keepPosterior = T, verbose=F, plot=F, cpus=5) {
  if ( coefs == 0 ) coefs = colnames(fit)
  if ( is.numeric(coefs) ) coefs = colnames(fit)[coefs]
  require(parallel)
  posts = posteriors(fit, coefs = coefs, quiet= !verbose, plot=plot, cpus=cpus, prior='empirical')
  ranks = postRank(posts, quiet = !verbose, cpus=cpus)

  if ( !keepPosterior ) fitSmall = fit
  fit$posterior = posts$data
  fit$prior = posts$prior
  fit$XRank = as.matrix(ranks$XRank)
  fit = bestGuess(fit, coefs = coefs, quiet = !verbose, cpus=cpus)
  
  if ( plot )  {
    n = length(fit$prior)
    m = round(sqrt(n))
    nr = floor(n/m)
    nc = ceiling(n/nr)
    layout(matrix(1:(nc*nr), nrow=nr, byrow=T))
    for ( i in 1:n )
      plotXRank(fit, coef=i, legend.cex=1)
    layout(1)
  }

  if ( !keepPosterior ) {
    fitSmall$prior = fit$prior
    fitSmall$XRank = fit$XRank
    fitSmall$best.guess = fit$best.guess
    fit = fitSmall
  }
  return(fit)
}


plotXRank = function(fit, top = 5, coef=0, verbose=F, legend.cex=1, ...) {
  if ( !('posterior' %in% names(fit)) )
    warning('Could not find posterior. Run XRank(fit).')
  if ( !('XRank' %in% names(fit)) )
    warning('Could not find XRank. Run XRank(fit).')
  if ( coef == 0 ) {
    if ( sum(fit$design[,1]) == length(fit$design[,1]) )
      coef = colnames(fit$coefficients)[2]
    else
      coef = colnames(fit$coefficients)[1]
  }
  if ( is.numeric(coef) )
    coef = names(fit$posterior)[coef]
  if ( verbose ) cat('Using contrast', coef, '\n')

  x = as.numeric(colnames(fit$posterior[[coef]]))
  binsize = fit$prior[[coef]]$breaks[2:(length(x)+1)] -
            fit$prior[[coef]]$breaks[1:length(x)]
  
  topGenes = order(fit$XRank[,coef])[1:top]
  topPosteriors = fit$posterior[[coef]][topGenes,]

  cols = c('red', 'orange', 'green', 'cyan', 'blue', 'violet', 'darkblue')
  if ( top > length(cols) ) cols = c(cols, rep('black', top-length(cols)))
  cols = mcri(cols)

  smooth = list()
  meanX = list()
  ymax = 0
  for ( i in 1:top ) {
    smooth[[i]] = spline(x, topPosteriors[i,]/binsize)
    meanX[[i]] = sum(smooth[[i]]$x*smooth[[i]]$y)/sum(smooth[[i]]$y)
    ymax = max(ymax, max(smooth[[i]]$y))
  }
  xmin = min(-0.5, 1.3*min(unlist(meanX)))
  xmax = max(0.5, 1.3*max(unlist(meanX)))

  
  plot(1,1, type='n', xlim=c(xmin, xmax), ylim=c(0, ymax),
       xlab = 'LFC', ylab = 'probability density',
       main = paste('top differential genes for', coef), ...)
  for ( i in 1:top ) {
    lines(smooth[[i]]$x, noneg(smooth[[i]]$y), col = cols[i], lwd=3)
    bg = fit$best.guess[topGenes[i],coef]
    ybg = interpolate(smooth[[i]]$x, smooth[[i]]$y, bg)
    segments(bg, -ymax/10, bg, ybg+ymax/10, col=cols[i], lwd=1)
  }

  if ( (xmax > abs(xmin) && xmin != -0.5) || xmax == 0.5 ) side = 'topright'
  else side = 'topleft'
  legend(side, paste(rownames(fit)[topGenes], ' (FDR ', signif(p.adjust(fit$p.value[,coef], method='fdr')[topGenes]),')', sep=''),
         col = cols[1:top], lwd=5, cex=legend.cex)
}

plotBayesPDFs = function(fit, coef = 0, ranks = 1) {
  if ( coef == 0 ) coef = ncol(fit$XRank)
  genes = rank(fit$XRank[,coef])[ranks]

  prior = fit$prior[[coef]]
  beta = fit$coefficients[,coef]
  w = fit$coefficients[,coef]/fit$t[,coef]
  df = fit$df.total[,]
  posterior = fit$posterior[[coef]]
  bg = fit$best.guess[,coef]
  
  doPlotBayes(prior=prior, beta=beta[genes[1]], w=w[genes[1]], df=df[genes[1]],
              posterior = posterior[genes[1]], bg=bg[genes[1]], add=F)
  if ( length(genes) > 1 ) {
    for ( gene in genes[-1] ) {
      doPlotBayes(prior=prior, beta=beta[gene], w=w[gene], df=df[gene],
                  posterior = posterior[gene], bg=bg[gene], add=T)
        
    }
  }
}

doPlotBayes = function(prior, beta, w, df, posterior, bg=NA, add=F) {
  
}


postWidth = function(fit) {
  
  postWidth = function(col) sapply(1:nrow(fit), function(gene)
    sqrt(sum(fit$prior[[col]]$mids^2*fit$posterior[[col]][gene,]) -
         sum(fit$prior[[col]]$mids*fit$posterior[[col]][gene,])^2))
  postWidths = do.call(cbind, lapply(1:ncol(fit), postWidth))
  colnames(postWidths) = colnames(fit)
  rownames(postWidths) = rownames(fit)
  
  fit$postWidth = postWidths
  return(fit)
}


topXRank = function(fit, coef=1, number=10) {
  fit = subsetFit(fit, rows=order(fit$XRank)[1:number], cols=coef)
  out = data.frame('gene' = rownames(fit), 'LFC'=as.numeric(fit$coefficient),
    'best.guess' = as.numeric(fit$best.guess), 'p.value' = as.numeric(fit$p.value),
    'FDR'=p.adjust(fit$p.value, method='fdr') ,'XRank' = as.numeric(fit$XRank), stringsAsFactors=F)
  invisible(out)
}
