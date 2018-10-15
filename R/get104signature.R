

#' returns the 104 signature from a superFreq q variant object
#'
#' @import MutationalPatterns
#' @import BSgenome
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
get104profile = function(q, Rdirectory, genome, somaticPcut = 0.5) {
  mx96 = get96signature(q=q, Rdirectory=Rdirectory, genome=genome, somaticPcut=somaticPcut)
  mxIndel = getIndelSignature(q, somaticPcut=somaticPcut)
  mx104 = rbind(mx96, mxIndel)
  return(mx104)
}


#get the 8 types of indel counts.
getIndelSignature = function(q, somaticPcut=0.5) {
  q = q[q$somaticP > somaticPcut,]
  ins1 = sum(grepl('^\\+', q$variant) & nchar(q$variant) == 2)
  ins2 = sum(grepl('^\\+', q$variant) & nchar(q$variant) %in% 3:4)
  ins4 = sum(grepl('^\\+', q$variant) & nchar(q$variant)  %in% 5:7)
  ins7p = sum(grepl('^\\+', q$variant) & nchar(q$variant) > 7)
  #the gsub throws a warning and returns NA when q$variant is not on '-[0-9]' format
  #those will turn FALSE from & grepl('^-', q$variant) though, so suppressing warnings here.
  del1 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) == 1 & grepl('^-', q$variant)))
  del2 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) %in% 2:3 & grepl('^-', q$variant)))
  del4 = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) %in% 4:6 & grepl('^-', q$variant)))
  del7p = suppressWarnings(sum(as.numeric(gsub('-', '', q$variant)) > 6 & grepl('^-', q$variant)))
  
  indelMx = matrix(c(ins1, ins2, ins4, ins7p, del1, del2, del4, del7p), ncol=1,
    dimnames=list(c('ins1', 'ins2', 'ins4', 'ins7p', 'del1', 'del2', 'del4', 'del7p'), 'name'))
  
  return(indelMx)
}


#get the 96 signature through mutationalPatterns
get96signature = function(q, Rdirectory, genome, somaticPcut=0.5) {

  if ( !dir.exists(Rdirectory) ) {
    stop('R directory doesnt exist for get96signature.')
  }
  
  ref_genome = genomeToMPgenome(genome)
  
  vcfFile = paste0(Rdirectory, '/temp_', runif(1, 0, 1e10), '.vcf')
  qSom = q[q$somaticP > somaticPcut & q$variant %in% c('A', 'T', 'C', 'G'),]

  noMuts = F
  if ( nrow(qSom) == 0 ) {
    noMuts = T
    qSom = q[1,]
  }
    
  superFreq:::writeToVCF(qSom, vcfFile, genome=genome)
  vcf = read_vcfs_as_granges(vcfFile, 'temp', genome = ref_genome)
    
  if ( !all(seqlevels(vcf[[1]]) %in% seqlevels(get(ref_genome))) )
    vcf = lapply(vcf, function(x) rename_chrom(x))
  
  mx = mut_matrix(vcf_list = vcf, ref_genome = ref_genome)
  if ( noMuts ) mx = mx*0

  unlink(vcfFile)
  
  return(mx)
}

#translates the superFreq assembly names to the BSgenome assembly names
genomeToMPgenome = function(genome) {
  if ( genome == 'hg19' ) return('BSgenome.Hsapiens.UCSC.hg19')
  if ( genome == 'hg38' ) return('BSgenome.Hsapiens.UCSC.hg38')
  if ( genome == 'mm10' ) return('BSgenome.Mmusculus.UCSC.mm10')
  if ( genome == 'mm9' ) return('BSgenome.Mmusculus.UCSC.mm9')
  warning('Couldnt find genome')
  return('')
}

#plot the 104 profile
plot_104_profile = function (mut_matrix, ymax = 0, legend=T, main=paste0(colnames(mut_matrix), ' (', sum(mut_matrix), ' mutations)'), ...) {
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x/pmax(1,sum(x)))
  context = c(MutationalPatterns:::TRIPLETS_96, 'ins1', 'ins2', 'ins4', 'ins7', 'del1', 'del2', 'del4', 'del7')
  substitution = c(rep(MutationalPatterns:::SUBSTITUTIONS, each = 16), rep('indel',8))
  substring(context[1:96], 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = reshape2::melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if ( ymax==0 ) ymax = max(df3$value, 0.01)
  cols = superFreq:::mcri(c('red', 'green', 'magenta', 'blue', 'cyan', 'orange', 'violet'))
  a=barplot(df3$value, space = c(0.4,ifelse(df3$substitution[1:103] == df3$substitution[2:104], 0.4, 2)),
          col=cols[as.numeric(as.factor(df3$substitution))], ylim=c(-ymax*0.2,ymax), main=main, ...)

  baseToCol = c('C'=superFreq:::mcri('red'), 'G'=superFreq:::mcri('orange'), 'T'=superFreq:::mcri('blue'), 'A'=superFreq:::mcri('cyan'))
  base1 = gsub('\\..+$', '', df3$context)[1:96]
  base2 = gsub('>.+$', '', df3$substitution)[1:96]
  base3 = gsub('^.+]', '', df3$context)[1:96]
  xs = a[1:96,1]
  rect(xleft=xs-0.5, ybottom=-0.04*ymax, xright=xs+0.5, ytop=-0.01*ymax, border=F, col=baseToCol[base1])
  rect(xleft=xs-0.5, ybottom=-0.07*ymax, xright=xs+0.5, ytop=-0.04*ymax, border=F, col=baseToCol[base2])
  rect(xleft=xs-0.5, ybottom=-0.10*ymax, xright=xs+0.5, ytop=-0.07*ymax, border=F, col=baseToCol[base3])

  text(a[97:104,1], -ymax*0.05, df3$context[97:104], cex=0.6, srt=90)
  
  for ( subs in unique(df3$substitution) ) {
    x = mean(a[df3$substitution == subs,])
    text(x, -ymax*0.15, subs, cex=1, col=cols[as.numeric(as.factor(unique(df3$substitution)))][subs==unique(df3$substitution)])
  }

  if ( legend )
    legend('topleft', names(baseToCol), pch=15, col=baseToCol, pt.cex=2, cex=0.8)
}


#plot the 104 profiles by sample
plot104profilesBySample = function(qs, Rdirectory, plotDirectory, genome, somaticPcut=0.5) {
  patternDirectory = paste0(plotDirectory, '/mutationalPatterns')
  superFreq:::ensureDirectoryExists(patternDirectory)
  pdf(paste0(patternDirectory, '/bySample.pdf'), width=10, height=6)
  mxList = lapply(names(qs), function(sample) {
    mx104 = get104profile(q=qs[[sample]], Rdirectory=Rdirectory, genome=genome, somaticPcut=somaticPcut)
    plot_104_profile(mx104, main=paste0(sample, ' (', sum(mx104), ' mutations)'))
    return(mx104)
  })
  dev.off()
  mx = do.call(cbind, mxList)
  colnames(mx) = names(qs)
  invisible(mx)
}

#plot the 104 profiles by clone
plot104profilesByClone = function(storyList, qs, Rdirectory, plotDirectory, genome) {
  patternDirectory = paste0(plotDirectory, '/mutationalPatterns')
  superFreq:::ensureDirectoryExists(patternDirectory)
  q = qs[[1]]
  pdf(paste0(patternDirectory, '/byClone.pdf'), width=10, height=6)
  mxList = lapply(names(storyList), function(clone) {
    mutations = storyList[[clone]]
    mutations = mutations[mutations %in% rownames(q)]
    if ( length(mutations) == 0 ) return()
    mx104 = superFreq:::get104profile(q=q[mutations,], Rdirectory=Rdirectory, genome=genome, somaticPcut=-1)
    colnames(mx104) = clone
    superFreq:::plot_104_profile(mx104, main=paste0('clone ', clone, ' (', sum(mx104), ' mutations)'))
    return(mx104)
  })
  dev.off()
  mx = do.call(cbind, mxList)
  invisible(mx)
}

#load the variants and stories from the R directory.
#calculate the profiles.
#plot the profiles by sample and clone.
#save the signatures to the Rdirectory.
loadAndPlotProfiles = function(Rdirectory, plotDirectory, genome) {
  storyFile = paste0(Rdirectory, '/stories.Rdata')
  variantFile = paste0(Rdirectory, '/allVariants.Rdata')
  if ( !all(file.exists(c(storyFile, variantFile))) ) {
    warning('Couldnt find data in ', Rdirectory, ' for mutational profiles. Could be wrong path or run hasnt finished.')
    return()
  }
  load(storyFile)
  stories = stories$stories
  load(variantFile)
  variants = allVariants$variants
  plotProfiles(Rdirectory=Rdirectory, plotDirectory=plotDirectory, stories=stories,
               variants=variants, genome=genome)
}

plotProfiles = function(Rdirectory, plotDirectory, stories, variants, genome, forceRedo=F) {
  catLog('Mutaional signatures of somatic mutations...')
  saveFile = paste0(Rdirectory, '/profiles.Rdata')
  if ( !file.exists(saveFile) | forceRedo ) {
    storyList = stories[[1]]$consistentClusters$storyList
    qs = variants$variants
    catLog('by sample..')
    sampleMx = plot104profilesBySample(qs=qs, Rdirectory=Rdirectory, plotDirectory=plotDirectory, genome=genome, somaticPcut=0.5)
    catLog('by clone..')
    cloneMx = plot104profilesByClone(storyList=storyList, qs=qs, Rdirectory=Rdirectory, plotDirectory=plotDirectory, genome=genome)
    profiles = list(sampleMx=sampleMx, cloneMx=cloneMx)
    catLog('saving results..')
    save(profiles, file=saveFile)
    catLog('done.\n')
  }
  else catLog('save file already there, moving along.\n')
}
