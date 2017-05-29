# superFreq
SuperFreq is an R package that analyses cancer exomes.

![](https://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/bfafbe20134792cfd1803e3381849e0a54aa5b9d/images/river.png)

# What does it do?
SuperFreq analyses and filters somatic SNVs and short indels, calls copy numbers and tracks clones over multiple samples from the same individual. It identifies the copy number alterations and point mutations in each clone, and highlights potentially causing mutations through variant annotation and COSMIC.

![](https://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/1b2575a930b5bf40b74d0a34ec4b708ebba89e1c/images/all.png)


# How do I run it?
Start R

```
library(devtools)
install_github('ChristofferFlensburg/superFreq')

library(superFreq)
?superFreq
```

A typical analysis first sets the parameters and then calls the `superFreq()` function:

```
library(superFreq)

cpus=12

#this file needs to be created. See ?superFreq
metaDataFile = 'metaData.tsv'

#a bed file with the capture regions of the exome. Or just the exones if capture regions are not available.
captureRegionsFile = '~/resources/captureRegions/myCaptureInThisBatch.bed'

#This directory needs to be created and set up. See ?superFreq
normalDirectory = '~/resources/superFreq/referenceNormals/myCaptureInThisBatch'

#The reference fasta and name. hg19, hg38 and mm10 available.
reference = '~/resources/reference/hg19/hg19.fa'
genome = 'hg19'

#the dbSNP and cosmic directory. This will be created and downloaded if not existing.
dbSNPdirectory = normalizePath('~/resources/superFreq/dbSNP')
cosmicDirectory='~/resources/superFreq/COSMIC'

#The directory where the log file and saved .Rdata is stored. Will be created.
Rdirectory = 'R'
#The directory where all the plots and tables from the analysis go. Will be created.
plotDirectory = 'plots'

#superFreq reuses saved data if available. This setting can force it to redo part of the analysis.
#default forceRedoNothing() means that saved information is used whenever available.
forceRedo = forceRedoNothing()

#a measure on how much large-scale biases are expected in the coverage.
#this controls the sensitivity vs accuracy of the coverage part of copy number calls.
systematicVariance=0.02
#a measure on how much biases (such as PCR duplication) is expected in the VAFs.
#this controls the sensitivity vs accuracy of the heterozygous SNP part of copy number calls.
maxCov=150

#The format of the quality scores of the base calls. Almost always 33 these days.
BQoffset = 33

#The mode. 'DNA' is for exomes, while 'RNA' has some minor changes when running on RNA.
mode = 'DNA'

#This setting runs each individual separately (as indicated in the metadata).
#will create subdirectories in the plotDirectory and Rdirectory.
#This is suggested whenever there is more than one individual in the batch.
splitRun = T

data =
    superFreq(metaDataFile, captureRegions=captureRegionsFile, normalDirectory=normalDirectory,
              normalCoverageDirectory=normalCoverageDirectory,
              Rdirectory=Rdirectory, plotDirectory=plotDirectory, reference=reference, genome=genome,
              BQoffset=BQoffset, cpus=cpus, forceRedo=forceRedo, systematicVariance=systematicVariance,
              maxCov=maxCov, mode=mode, splitRun=splitRun)
```


More information is in the manual. 

![](https://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/bfafbe20134792cfd1803e3381849e0a54aa5b9d/images/multisample.png)

# What is the input?
You need the aligned bam files of the exomes, and a preliminary (liberal) variant calling through for example varScan, mutect, multiSNV or any other similar software. superFreq is not sensitive to false positives in these VCFs.

SuperFreq also requires a set of (at least 2, 5+ is better) reference normal samples from the same sequencing technology and exome capture.
Preferably sequenced in the same lab. These samples do not have to be related to the analysed cancer exomes.

You also need some meta data:
- a .bed file with the capture regions
- the fasta file you aligned to
- a tab separated file with information about the samples

![](https://gitlab.wehi.edu.au/flensburg.c/superFreq/raw/f7e93765a81a7dd349d5095a8aba16fecf96b43c/images/TCGA.A3.3320.PrimaryTumor.WXS.fc2.png)

# What is the output?
Plots (some shown here), tables, and R objects for downstream analysis. Analysis results as well as diagnostic quality control. Some examples:
- Scatter plots of variant frequencies between pairs of samples
- heatmaps of somatic variant frequencies
- CNA plots for each samples
- river plots showing clonal evolution
- summary plots showing SNVs and CNAs over all samples
- top tables of annotated somatic variants.

#dependencies
- R.
- a bunch of R packages.
- VEP version 76 or later.
