# superFreq
SuperFreq is an R package that analyses cancer exomes.

# What does it do?
SuperFreq analyses and filters somatic SNVs and short indels, calls copy numbers and tracks clones over multiple samples from the same individual. It identifies the copy number alterations and point mutations in each clone, and highlights potentially causing mutations through variant annotation and COSMIC.

# How do I run it?
Start R

```
install.packages('devtools')
library(devtools)
install_github('ChristofferFlensburg/superFreq')

?superFreq
```

More information is in the manual. 

# What is the input?
You need the aligned bam files of the exomes, and a preliminary (liberal) variant calling through for example varScan, mutect, multiSNV or any other similar software. superFreq is not sensitive to false positives in these VCFs.

SuperFreq also requires a set of (at least 2, 5+ is better) reference normal samples from the same sequencing technology and exome capture.
Preferably sequenced in the same lab. These samples do not have to be related to the analysed cancer exomes.

You also need some meta data:
- a .bed file with the capture regions
- the fasta file you aligned to
- a tab separated file with information about the samples

# What is the output?
Plots, tables, and R objects for downstream analysis. Analysis results as well as diagnostic quality control. Some examples:
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
