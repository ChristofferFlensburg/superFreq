
#Compares all pairs of samples from the same individual and outputs signidicantly different frequencies between the two samples.
#These are essentially the red points in the frequency scatter plots.
outputNewVariants = function(variants, pairs, genome, directory, cpus=1, forceRedo=F) {
  outfile = paste0(directory, '/newVariants.xls')
  if ( (!file.exists(outfile) | forceRedo) & length(pairs) > 0 ) {
    news = list()
    for ( pair in pairs ) {
      name1 = substring(gsub('\\.', '', paste0(pair[1], ' to ', pair[2])), 1, 31)
      name2 = substring(gsub('\\.', '', paste0(pair[2], ' to ', pair[1])), 1, 31)
      catLog('Looking for new cancer variants in ', name1, '\n')
      news[[name1]] = newVariants(variants$variants[[pair[1]]], variants$variants[[pair[2]]], variants$SNPs, genome, cpus=cpus)
      catLog('Looking for new cancer variants in ', name2, '\n')
      news[[name2]] = newVariants(variants$variants[[pair[2]]], variants$variants[[pair[1]]], variants$SNPs, genome, cpus=cpus)
    }
    WriteXLS('news', outfile)
  }
}

#helper function that isolates the significantly different variants
newVariants = function(q1, q2, SNPs, genome='hg19', cpus=1, ps=NA) {
  common = intersect(rownames(q1), rownames(q2))
  common = common[!is.na(common)]
  single = rownames(q2)[!(rownames(q2) %in% rownames(q1))]
  if ( length(common) == 0 ) return(data.frame())
    
  q1 = q1[common,]
  q2 = q2[common,]
  freq1 = q1$var/q1$cov
  freq1[is.na(freq1)] = -0.02
  freq2 = q2$var/q2$cov
  freq2[is.na(freq2)] = -0.02
  
  if ( is.na(ps[1]) ) {
    if ( cpus==1 )   
      ps = sapply(1:length(freq1), function(i) fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value)
    else {
      ps = unlist(mclapply(1:length(freq1), function(i) fisher.test(matrix(c(q1$ref[i], q1$var[i], q2$ref[i], q2$var[i]), nrow=2))$p.value, mc.cores=cpus))      
    }
  }
  
  low = (freq1 < 0.2 | (pBinom(q1$cov, q1$var, 0.3) < 0.01 & freq1 < 0.3)) & !(freq1 == 0 & freq2 == 0)
  fdr = p.adjust(ps, method='fdr')
  fdr[low] = p.adjust(ps[low], method='fdr')
  use = which(freq1 < freq2 & low)

  cut = 1/sum(low)^0.75
  Nprint = sum(ps[use] <= cut) + 10
  toReturn = use[order(ps[use])][1:Nprint]
  severity =
    if ( 'severity' %in% names(q1) & 'severity' %in% names(q2) )
      pmin(q1$severity[toReturn], q2$severity[toReturn])
    else
      rep('na', nrow(q1[toReturn,]))
  effect =
    if ( 'type' %in% names(q1) & 'type' %in% names(q2) )
      ifelse(q1$severity[toReturn] == 100, q2$type[toReturn], q1$type[toReturn])
    else
      rep('na', nrow(q1[toReturn,]))
  ret = data.frame(
    chr=xToChr(q1$x[toReturn], genome),
    start=xToPos(q1$x[toReturn], genome),
    end=xToPos(q1$x[toReturn], genome),
    reference=q1$reference[toReturn],
    variant=q1$variant[toReturn],
    inGene=SNPs[as.character(q1$x[toReturn]),]$inGene,
    severity=severity,
    effect= effect,
    f1=q1$var[toReturn]/q1$cov[toReturn],
    f2=q2$var[toReturn]/q2$cov[toReturn],
    cov1=q1$cov[toReturn],
    ref1=q1$ref[toReturn],
    var1=q1$var[toReturn],
    cov2=q2$cov[toReturn],
    ref2=q2$ref[toReturn],
    var2=q2$var[toReturn],
    pDiff=ps[toReturn],
    fdr=fdr[toReturn],
    flag1=q1$flag[toReturn],
    flag2=q2$flag[toReturn],
    pbq1=q1$pbq[toReturn],
    pbq2=q2$pbq[toReturn],
    pmq1=q1$pmq[toReturn],
    pmq2=q2$pmq[toReturn],
    psr1=q1$psr[toReturn],
    psr2=q2$psr[toReturn],
    dbSNP=ifelse(q1$db[toReturn], 'YES', ''),
    germlineLike=ifelse(is.na(q1$germline[toReturn]), 'na', ifelse(q1$germline[toReturn], 'YES', '')),
    row.names=rownames(q1)[toReturn])

  if ( all(moreVEPnames(genome=genome) %in% names(q2)) ) {
    catLog('Adding more VEP info..')
    ret = cbind(ret, q2[toReturn,moreVEPnames(genome=genome)])
    catLog('done.\n')
  }

  
  if ( 'severity' %in% names(q1) ) ord = order(q1$severity[toReturn] + 10*q1$germline[toReturn])
  else ord = order(10*q1$germline[toReturn])
  ret = ret[ord,]
  return(ret)
}




