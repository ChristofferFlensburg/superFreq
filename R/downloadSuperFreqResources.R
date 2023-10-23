

downloadSuperFreqDbSNP = function(directory='superFreqDbSNP', genome='hg19') {
  if ( genome == 'hg19' ) {
    if ( !dir.exists(directory) ) {
      cat('Creating superFreq dbSNP annotation directory ', directory, '.\n')
      dir.create(directory)
    }

    #hg19 dbSNP
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/dbAFnew.Rdata')
	figshareUrl = 'https://figshare.com/ndownloader/files/31278718'
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)

    #hg19 ExAC
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/exac.Rdata')
	figshareUrl = 'https://figshare.com/ndownloader/files/31278721'
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)
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
	figshareUrl = 'https://figshare.com/ndownloader/files/31278724'
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)

    #hg19 ExAC
    url = paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/dbSNP/hg38/exac.Rdata')
	figshareUrl = 'https://figshare.com/ndownloader/files/31278727'
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)
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
	figshareUrl = 'https://figshare.com/ndownloader/files/31278730'
    destFile = paste0(directory, '/', basename(url))
    if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)
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
    filenames = basename(urls)
	figshareUrls = c('https://figshare.com/ndownloader/files/31278697',
					 'https://figshare.com/ndownloader/files/31278712',
					 'https://figshare.com/ndownloader/files/31278706')
    if ( genome == 'hg38' ) {
        urls =
            paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/',
                   c('COSMIC/allCosmicCounts_hg38.Rdata', 'COSMIC/cosmicCounts.Rdata', 'ClinVar/clinvar38.Rdata'))
       filenames = basename(urls)
	   figshareUrls = c('https://figshare.com/ndownloader/files/31278700',
						'https://figshare.com/ndownloader/files/31278712',
						'https://figshare.com/ndownloader/files/31278709')
 
    }
    if ( genome == 'mm10' ) {
      urls =
        c('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/COSMIC/CCGD_export.csv',
          'https://figshare.com/ndownloader/files/36269079')
      filenames = c('CCGD_export.csv', 'MGIBatchReport_20220629_032731_filtered.txt')
	  figshareUrls = c('https://figshare.com/ndownloader/files/31278703',
					   'https://figshare.com/ndownloader/files/36269079')
    }
    for ( i in 1:length(figshareUrls) ) {
      url = figshareUrls[i]
      destFile = paste0(directory, '/', filenames[i])
      if ( !file.exists(destFile) ) superFreq:::downloadFile(url, destfile = destFile)
    }
  
}

#change how this is called to download to superFreq resources rather than to the R directory.
downloadCaptureRegions = function(directory, genome='hg19', mode='exome', binSize=10000) {
    url = ''
    if ( genome == 'hg19' & mode %in% c('exome', 'RNA') ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg19exons.bed'
		figshareUrls = c('https://figshare.com/ndownloader/files/31278643',
						 'https://figshare.com/ndownloader/files/31278646',
						 'https://figshare.com/ndownloader/files/31278649')
    }
    if ( genome == 'hg38' & mode %in% c('exome', 'RNA') ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg38exons.bed'        
		figshareUrls = c('https://figshare.com/ndownloader/files/31278661',
						 'https://figshare.com/ndownloader/files/31278664',
						 'https://figshare.com/ndownloader/files/31278667')
    }
    if ( genome == 'mm10' & mode %in% c('exome', 'RNA') ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/mm10exons.bed'
		figshareUrls = c('https://figshare.com/ndownloader/files/31278679',
						 'https://figshare.com/ndownloader/files/31278682',
						 'https://figshare.com/ndownloader/files/31278685')
    }

    if ( genome == 'hg19' & mode == 'genome' & binSize == 10000 ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg19genome.10k.bed'
		figshareUrls = c('https://figshare.com/ndownloader/files/31278652',
						 'https://figshare.com/ndownloader/files/31278655',
						 'https://figshare.com/ndownloader/files/31278658')
	}
    if ( genome == 'hg38' & mode == 'genome' & binSize == 10000 ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/hg38genome.10k.bed'
		figshareUrls = c('https://figshare.com/ndownloader/files/31278670',
						 'https://figshare.com/ndownloader/files/31278673',
						 'https://figshare.com/ndownloader/files/31278676')
	}
    if ( genome == 'mm10' & mode == 'genome' & binSize == 10000 ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/captureRegions/mm10genome.10k.bed'
		figshareUrls = c('https://figshare.com/ndownloader/files/31278688',
						 'https://figshare.com/ndownloader/files/31278691',
						 'https://figshare.com/ndownloader/files/31278694')
	}

    if ( url == '' ) {
        catLog('\nERROR: Couldnt find any capture regions for genome ', genome, ' in mode ', mode, ' with binSize ', binSize, '\n\n')
        stop('Couldnt find any capture regions for genome ', genome, ' in mode ', mode, ' with binSize ', binSize)
    }
    
    urls = c(url, gsub('.bed$', '.named.bed', url), gsub('.bed$', '.named.dn.bed', url))
    for ( i in 1:3 ) {
      url = urls[i]
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) downloadFile(figshareUrls[i], destfile = destFile)
    }

    return(paste0(directory, '/', basename(urls[1])))
}


downloadSuperFreqAnnotation = function(directory='superFreqAnnotation', genome='hg19') {
    if ( !dir.exists(directory) ) {
      cat('Downloading superFreq annotation to', directory, '.\n')
      dir.create(directory)
    }

    url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/dumpHg19.Rdata'
    figshareUrl = 'https://figshare.com/ndownloader/files/31278631'
    if ( genome == 'hg38' ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/dumpHg38.Rdata'
		figshareUrl = 'https://figshare.com/ndownloader/files/31278634'        
    }
    if ( genome == 'mm10' ) {
        url = 'http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/annotation/dumpMm10.Rdata'
		figshareUrl = 'https://figshare.com/ndownloader/files/31278637'        
    }

	destFile = paste0(directory, '/', basename(url))
	if ( !file.exists(destFile) ) downloadFile(figshareUrl, destfile = destFile)
  
}


downloadSuperFreqSignatures = function(directory='superFreqSignatures') {
    if ( !dir.exists(directory) ) {
      cat('Downloading superFreq signatures to', directory, '.\n')
      dir.create(directory)
    }

    urls =
        paste0('http://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/master/',
               c('signatures/extendedPatternMx.Rdata', 'signatures/cosmicSignatures_v2_104.Rdata',
                 'signatures/mutagenSignatures.Rdata'))
	figshareUrls = c('https://figshare.com/ndownloader/files/31278736',
					 'https://figshare.com/ndownloader/files/31278733',
					 'https://figshare.com/ndownloader/files/31278739')

    for ( i in 1:3 ) {
      url = urls[i]
      destFile = paste0(directory, '/', basename(url))
      if ( !file.exists(destFile) ) downloadFile(figshareUrls[i], destfile = destFile)
    }

}

downloadFile = function(url, destfile) {
	system(paste0('wget --no-check-certificate -O ', destfile, ' ', url))
}




