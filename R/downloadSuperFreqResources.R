

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
      paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/',
             c('COSMIC/allCosmicCounts_hg19.Rdata', 'COSMIC/cosmicCounts.Rdata', 'ClinVar/clinvar19.Rdata'))
    if ( genome == 'hg38' )
        urls =
            paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/',
                   c('COSMIC/allCosmicCounts_hg38.Rdata', 'COSMIC/cosmicCounts.Rdata', 'ClinVar/clinvar38.Rdata'))
    if ( genome == 'mm10' )
        urls = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/CCGD_export.csv'
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }
  
}

downloadCaptureRegions = function(directory, genome='hg19', mode='exome', binSize=10000) {
    url = ''
    if ( genome == 'hg19' & mode %in% c('exome', 'RNA') )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg19exons.bed'
    if ( genome == 'hg38' & mode %in% c('exome', 'RNA') )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg38exons.bed'
    if ( genome == 'mm10' & mode %in% c('exome', 'RNA') )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/mm10exons.bed'

    if ( genome == 'hg19' & mode == 'genome' & binSize == 10000 )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg19genome.10k.bed'
    if ( genome == 'hg38' & mode == 'genome' & binSize == 10000 )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg38genome.10k.bed'
    if ( genome == 'mm10' & mode == 'genome' & binSize == 10000 )
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/mm10genome.10k.bed'

    if ( url == '' ) {
        catLog('\nERROR: Couldnt find any capture regions for genome ', genome, ' in mode ', mode, ' with binSize ', binSize, '\n\n')
        stop('Couldnt find any capture regions for genome ', genome, ' in mode ', mode, ' with binSize ', binSize)
    }
    
    urls = c(url, gsub('.bed$', '.named.bed', url), gsub('.bed$', '.named.dn.bed', url))
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }

    return(paste0(directory, '/', basename(urls[1])))
}


downloadSuperFreqAnnotation = function(directory='superFreqAnnotation', genome='hg19') {
    if ( !dir.exists(directory) ) {
      cat('Downloading superFreq annotation to', directory, '.\n')
      dir.create(directory)
    }

    urls =
        paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/',
               c('dumpHg19.Rdata'))
    if ( genome == 'hg38' )
        urls =
            paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/',
                   c('dumpHg38.Rdata'))
    if ( genome == 'mm10' )
        urls =
            paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/',
                   c('dumpMm10.Rdata'))
    for ( url in urls ) {
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) download.file(url, destfile = destFile, method="wget")
    }
  
}


