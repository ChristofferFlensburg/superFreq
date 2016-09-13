

downloadSuperFreqDbSNP = function(directory='superFreqDbSNP', genome='hg19') {
  if ( genome == 'hg19' ) {
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation directory ', directory, '.\n')
      dir.create(directory)
    }

    #hg19 dbSNP
    chrs = c(1:22, 'X', 'Y', 'MT')
    urls = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/ds_flat_ch', chrs, '.Rdata')
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }

    #hg19 ExAC
    chrs = c(1:22, 'X', 'Y')
    urls = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/', chrs, '.Rdata')
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }
  }
  else stop('Genome ', genome, ' not supported for dbSNP in superFreq, sorry.')
}


downloadSuperFreqCOSMIC = function(directory='superFreqCOSMIC', genome='hg19') {
  if ( genome == 'hg19' ) {
    if ( !dir.exists(directory) ) {
      cat('Downloading superFreq COSMIC annotation to', directory, '.\n')
      dir.create(directory)
    }
    
    urls =
      paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/',
             c('allCosmicCounts', 'cosmicCounts'), '.Rdata')
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }
  }
  else stop('Genome ', genome, ' not supported for dbSNP in COSMIC, sorry.')
  
}
