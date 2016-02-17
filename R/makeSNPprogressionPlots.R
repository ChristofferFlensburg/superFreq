

#prints heatmaps and line plots of the frequency of the SNVs in the samples of the same individual
#also outputs the results to excel files.
makeSNPprogressionPlots = function(variants, timeSeries, normals, plotDirectory, cpus=1, forceRedo=F) {
  msDirectory = paste0(plotDirectory, '/multiSample')
  if ( length(timeSeries) == 0 ) return()
  if ( !file.exists(msDirectory) ) dir.create(msDirectory)
  for ( i in 1:length(timeSeries) ) {
    ts = timeSeries[[i]]
    outfile = paste0(msDirectory, '/', names(timeSeries)[i], '.pdf')
    excelFileDB = paste0(msDirectory, '/', names(timeSeries)[i], '_DB.xls')
    excelFileNotDB = paste0(msDirectory, '/', names(timeSeries)[i], '_NotDB.xls')
    excelFileAllChanging = paste0(msDirectory, '/', names(timeSeries)[i], '_allSevere.xls')
    if ( !file.exists(outfile) | forceRedo ) {
      catLog('Plotting SNP progression to ', outfile, '...\n', sep='')
      pdf(outfile, width = 15, height=10)
      qualityProgression(variants$variants[ts], variants$SNPs, normals[ts], nondb=F, excelFile=excelFileDB, main='significantly changing germline SNPs')
      qualityProgression(variants$variants[ts], variants$SNPs, normals[ts], db=F, excelFile=excelFileNotDB, main='significantly changing somatic SNVs')
      qualityProgression(variants$variants[ts], variants$SNPs, normals[ts], db=F, excelFile=excelFileAllChanging, main='all protein changing somatic SNVs', filterConstant=F)
      dev.off()
      catLog('done.\n')
    }
  }
}


#helper function that does all the work for the frequency progression plots.
qualityProgression = function(qs, SNPs, normal, db=T, nondb=T, excelFile='', main='', colMode='default', nCol=200, filterConstant=T, linePlot=T) {
  if ( !filterConstant  && !(('severity' %in% names(qs[[1]])) && !any(is.na(qs[[1]]$severity))) ) {
    catLog('Skipping unfiltered snp progression without VEP information.\n')
    return()
  }
  catLog('Finding common variants..')
  if ( !('somaticP' %in% names(qs[[1]])) ) {
    warning('no somaticQ in qs of quality progression plot. not plotting.')
    return()
  }
  if ( any(sapply(qs, function(q) any(is.na(q$somaticP)))) ) {
    warning('found NA somaticQ in qs of quality progression plot. not plotting.')
    return()
  }
  if ( !nondb ) use = lapply(qs, function(q) q$db & q$somaticP <= 0.1 & q$var > 0)
  if ( !db ) use = lapply(qs[!normal], function(q) q$somaticP > 0.1 & q$var > q$cov*0.1 & q$flag=='')
  use = apply(do.call(cbind, use), 1, any)
  if ( sum(use) == 0 ) return()
  qs = lapply(qs, function(q) q[use,] )
  if ( nondb & any(normal) ) {
    noisyNormals = unique(unlist(lapply(qs[normal], function(q) which(!q$db & q$cov > 0 & q$var/q$cov > 0.05))))
    if ( length(noisyNormals) > 0 )
      qs = lapply(qs, function(q) q[-noisyNormals,])
  }
  if ( nrow(qs[[1]]) == 0 ) return()
  catLog(nrow(qs[[1]]), '. Cleaning..', sep='')
  var = do.call(cbind, lapply(qs, function(q) q$var))
  keep = rowSums(var) > 0
  qs = lapply(qs, function(q) q[keep,])
  fs = do.call(cbind, lapply(qs, function(q) q$var/q$cov))
  colnames(fs) = names(qs)
  cov = do.call(cbind, lapply(qs, function(q) q$cov))
  var = do.call(cbind, lapply(qs, function(q) q$var))
  f = rowSums(var)/rowSums(cov)
  fs[is.na(fs)] = -0.02

  if ( filterConstant ) {
    catLog('p-values..')
    ps = matrix(pBinom(as.integer(cov), as.integer(var), rep(f, ncol(cov))), ncol=ncol(cov))
    p = apply(ps, 1, fisherTest)[2,]
  }
  else {
    severity = do.call(cbind, lapply(qs, function(q) q$severity))
    severity = apply(severity, 1, min)
    p = ifelse(severity <= 11, 0, 1)
  }

  SNPs = SNPs[SNPs$x %in% qs[[1]]$x,]
  gene = paste0(gsub('.+:', '', SNPs[as.character(qs[[1]]$x),]$inGene), ' (', SNPs[as.character(qs[[1]]$x),]$chr, ')')
  if ( db ) gene = SNPs[as.character(qs[[1]]$x),]$chr
  
  catLog('colours..')
  dof = max(20, sum(rowSums(fs) > 0 & rowSums(fs) < ncol(fs)))
  importance = pmin(1, noneg(-log10(p)/log10(dof) - 0.75))
  weight = pmin(1, pmax(importance, sqrt(rowMeans(cov/300))))
  doColour = importance > 0.5
  recurringGenes = gene[doColour][duplicated(gene[doColour])]
  isRecurringGene = gene %in% recurringGenes & doColour
  hue = rep(0, length(doColour))
  hue[isRecurringGene] = as.integer(as.factor(gene[isRecurringGene]))
  hue = hue/(max(hue)+1)
  col = D3colours(weight[doColour], importance[doColour], hue[doColour])

  if ( sum(doColour) > 1 ) {
    catLog('plotting heatmap..')
    nRecurringGenes = length(unique(hue[isRecurringGene]))
    rGcol = D3colours(rep(1, nRecurringGenes), rep(1, nRecurringGenes), unique(hue[isRecurringGene]))
    rG = as.character(unique(gene[isRecurringGene]))
    names(rGcol) = rG
    RSC = ifelse(gene[doColour] %in% rG, rGcol[gene[doColour]], 'grey')
    fs[fs==-0.02] = NA
    Rowv = NULL
    heatCol = colourGradient(cols=mcri(c('black', 'grey', 'cyan', 'blue', 'red')),
      anchors=c(0, 0.02, 0.2, 0.5, 1), steps=nCol)
    if ( colMode == 'sunset' ) {
        heatCol = colourGradient(cols=mcri(c('black', 'blue', 'cyan', 'orange')),
          anchors=c(0, 0.1, 0.5, 1), steps=nCol)
      if ( any(is.na(fs[doColour,,drop=F])) ) {
        nCol=1000
        fs[which(is.na(fs),arr.ind=T)] = 1.002*min(fs[doColour,,drop=F], na.rm=T) - 0.002*max(fs[doColour,,drop=F], na.rm=T)
        heatCol = colourGradient(cols=mcri(c('grey15', 'black', 'blue', 'cyan', 'orange')),
          anchors=c(0, 0.001, 0.1, 0.5, 1), steps=nCol)
      }
    }
    if ( nrow(fs[doColour,,drop=F]) <= 1000 ) {
      clusterOrder =
        try(makeHeatmap(fs[doColour,,drop=F], cexCol=1, labRow=gene[doColour], Rowv=Rowv, nCol=nCol, col=heatCol,
                        RowSideColors = RSC, margins=c(8,15), main=main, label='frequency'))
      if ( class(clusterOrder) == 'try-error' ) {
        Rowv = NA
        catLog('too many variants, not clustering...')
        clusterOrder =
          makeHeatmap(fs[doColour,,drop=F], cexCol=1, labRow=gene[doColour], Rowv=Rowv, nCol=nCol, col=heatCol,
                      RowSideColors = RSC, margins=c(8,15), main=main, label='frequency')
      }
    }
    else {
      Rowv = NA
      catLog('too many variants, not clustering...')
      clusterOrder =
        makeHeatmap(fs[doColour,,drop=F], cexCol=1, labRow=gene[doColour], Rowv=Rowv, nCol=nCol, col=heatCol,
                    RowSideColors = RSC, margins=c(8,15), main=main, label='frequency')
    }
    fs[is.na(fs)] = -0.02
    if ( length(rG) > 0 ) {
      legend('right', rG, col = rGcol, lwd=10, bg='white')
    }
    fs = fs[,clusterOrder[[2]]]
    N = ncol(fs)

    if ( linePlot ) {
      catLog('frequency progression..')
      plot(0, type='n', xlim=c(1, length(qs)*1.2), ylim=c(0,1), main=main)
      segments(col(fs)[doColour, 1:(N-1)], fs[doColour,1:(N-1)],
               col(fs)[doColour, 2:N], fs[doColour, 2:N],
               lwd=(weight[doColour]+importance[doColour]), col=RSC)
      text(1:N, 1.02, colnames(fs), cex=0.7)
      if ( length(rG) > 0 )
        legend('right', rG, col = rGcol, lwd=10, bg='white')
    }
    catLog('done!\n')
    
    if ( excelFile != '' ) {
      catLog('Output plotted variants to', excelFile, '...')    
      multiSampleData = data.frame(
        'gene'=gsub('.+:', '', SNPs[as.character(qs[[1]]$x[doColour]),]$inGene),
        'chr'=SNPs[as.character(qs[[1]]$x[doColour]),]$chr,
        'start'=xToPos(qs[[1]]$x[doColour]),
        'end'=xToPos(qs[[1]]$x[doColour]),
        'reference'=qs[[1]]$reference[doColour],
        'variant'=qs[[1]]$variant[doColour],
        'f'=fs[doColour,],
        'cov'=cov[doColour,clusterOrder[[2]]],
        'var'=var[doColour,clusterOrder[[2]]])
      
      if ( "severity" %in% names(qs[[1]]) ) {
        severity = do.call(cbind, lapply(qs, function(q) q$severity[doColour]))
        effect = do.call(cbind, lapply(qs, function(q) q$type[doColour]))
        colToUse = apply(severity, 1, function(qSeverities) which(min(qSeverities) == qSeverities)[1])
        severity = severity[cbind(1:nrow(severity), colToUse)]
        effect = effect[cbind(1:nrow(effect), colToUse)]
        multiSampleData = cbind(multiSampleData[,1:6], 'severity'=severity, 'effect'=effect,
          multiSampleData[,7:ncol(multiSampleData)])
      }
      
      write.csv(multiSampleData, gsub('.xls$', '.csv', excelFile))
      
      if ( nrow(multiSampleData) <= 65000 )
        WriteXLS('multiSampleData', excelFile)
      else {
        multiSampleData = multiSampleData[order(multiSampleData$severity)[1:65000],]
        WriteXLS('multiSampleData', excelFile)
      }
      
      
      catLog('done!\n')
    }
  }
  else catLog('skipping plots (no important variants).\n')
  
}

#helper function that calculates p-values for a sum of two binomials symmetric around 0.5
pTwoBinom = function(cov, var, f) {
  use = cov > 0 & var >= 0
  p = rep(1, length(cov))
  if ( length(f) == 1 ) f = rep(f, length(cov))
  cov = cov[use]
  var = var[use]
  f = f[use]
  if ( length(f) != length(cov) | length(var) != length(cov) ) cat('Length of f must match length of cov, or be 1.\n')
  f = pmin(f, refBiasMirror(f))
  fM = refBiasMirror(f)
  var = mirrorDown(var, cov)

  midPoints = twoBinom(round(cov*refBias(0.5)), cov, f, fM)
  density =twoBinom(var, cov, f, fM)
  outside = density < midPoints & var < cov*refBias(0.5)
  p = rep(1, length(var))
  if ( any(outside) ) {
    vo = var[outside]
    co = cov[outside]
    fo = f[outside]
    foM = fM[outside]
    p[outside] = pbinom(vo, co, fo)+pbinom(vo, co, foM) - twoBinom(vo, co, fo, foM)
  }
  if ( any(!outside) ) {
    vi = var[!outside]
    ci = cov[!outside]
    fi = f[!outside]
    viMirror = findMirror(vi, ci, fi)
    viLow = pmin(vi, viMirror)
    viHigh = pmax(vi, viMirror)
    p[!outside] = (pbinom(viLow, ci, fi) + pbinom(viLow, ci, 1-fi)) - (dbinom(viLow, ci, fi)+dbinom(viLow, ci, 1-fi))/2 +
      (pbinom(ci-viHigh, ci, fi) - pbinom(viHigh, ci, fi)) - (dbinom(ci-viHigh, ci, fi) - dbinom(viHigh, ci, fi))/2
  }
  p = pmin(1, p)    #rounding can give p-values just above 1. Move these down to 1.
  return(p)
}

#helper function, sum of two binomials
twoBinom = function(V, C, F1, F2) (dbinom(V, C, F1) + dbinom(V, C, F2))/2

#helper function, finding the other root of a sum of two binomials minus a constant.
findMirror = function(var, cov, f) {
  fM = refBiasMirror(f)
  #first guess for the mirror
  m = round(pmax(0, pmin(cov*refBias(0.5), 2*f*cov - var)))
  y = twoBinom(m, cov, f, fM)
  #target probability density to reach
  y0 = twoBinom(var, cov, f, fM)
  converged = rep(F, length(var))
  solutions = rep(0, length(var))
  while ( any(!converged) ) {
    tooLow = y < y0
    mnew = m - sign(y-y0)*sign(var-f*cov)
    ynew = twoBinom(mnew, cov, f, fM)
    solved = sign(y-y0) != sign(ynew-y0) | y == y0 

    #handle converged cases
    if ( any(solved) ) {
      solutions[which(!converged)[solved]] = ifelse(abs(y-y0)[solved] < abs(ynew-y0)[solved], m[solved], mnew[solved])
      converged[which(!converged)[solved]] = rep(T, sum(solved))
    }

    #update not converged cases, removing solved values
    m = mnew[!solved]
    y = ynew[!solved]
    y0 = y0[!solved]
    cov = cov[!solved]
    var = var[!solved]
    f = f[!solved]
    fM = fM[!solved]
  }
  return(solutions)
}

#helper function that picks a colour based on statistical weight, importance and hue.
#weight sets opaqueness, importance sets saturation, and hue sets the colour.
D3colours = function(weights, importances, hues) {
  return(apply(cbind(weights, importances, hues), 1, D3colour))
}

D3colour = function(D3) {
  weight = D3[1]
  importance = D3[2]
  hue = D3[3]
  x = 6*hue
  rgb = pmin(1, c(noneg(2-x)+noneg(x-4), noneg(x-noneg(2*x-4)), noneg(x-2-noneg(2*x-8))))
  rgb = 0.1*(1-importance) + rgb*importance
  return(rgb(rgb[1], rgb[2], rgb[3], weight))
}
