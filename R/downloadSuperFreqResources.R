

downloadSuperFreqDbSNP = function(directory='superFreqDbSNP', genome='hg19') {
  if ( genome == 'hg19' ) {
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation directory ', directory, '.\n')
      dir.create(directory)
    }

    #hg19 dbSNP
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/dbAFnew.Rdata')
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")

    #hg19 ExAC
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/exac.Rdata')
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
  }
  else if ( genome == 'hg38' ) {
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation directory ', directory, '.\n')
      dir.create(directory)
    }
    directory = paste0(directory, '/hg38')
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation subdirectory ', directory, '.\n')
      dir.create(directory)
    }

    #hg38 dbSNP
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/hg38/dbAFnew.Rdata')
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")

    #hg19 ExAC
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/hg38/exac.Rdata')
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
  }
  else if ( genome == 'mm10' ) {
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation directory ', directory, '.\n')
      dir.create(directory)
    }
    directory = paste0(directory, '/mm10')
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation subdirectory ', directory, '.\n')
      dir.create(directory)
    }

    #mm10 dbSNP
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/mm10/dbAFnew.Rdata')
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
  }
  else stop('Genome ', genome, ' not supported for dbSNP in superFreq, sorry.')
}


downloadSuperFreqCOSMIC = function(directory='superFreqCOSMIC', genome='hg19') {
    if ( !dir.exists(directory) ) {
      cat('Downloading superFreq COSMIC annotation to', directory, '.\n')
      dir.create(directory)
    }
    
    urls =
      paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/',
             c('allCosmicCounts.Rdata', 'cosmicCounts.Rdata', 'CCGD_export.csv'))
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }
  
}
