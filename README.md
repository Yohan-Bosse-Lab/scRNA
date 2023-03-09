## Overview  
 * This repository contains scripts that are used in the analysis of single-cell RNAseq data from lung adenocarcinoma patients
 * It contains an `Rmarkdown` script than can be run as such from the command line (or alternatively, `Knit` directly in `Rstudio`) :
 * Use at own risk (it's mainly at the Proof-Of-Concept stage).


## Installation
  * This is R-based, libraries requires are described in `.Rmd` files
```
library(dplyr) #data handling
library(Seurat) #scRNA data processing
library(patchwork) #data viz
library(sctransform) #data normalization
library(ggplot2) #data viz
library(DT) #data viz
library(scAnnotatR) #cell type annotation
```


## Basic usage
``` bash
#With some (default) parameters.
Rscript -e "rmarkdown::render('scRNA_seurat_lung_adenocarcino_vSCT.Rmd',params = list(nFeature_min=250, 
nCount_min=500, nCount_max=50000, percent_MT_max=3)"
```

## Further information
  * sebastien.renaut.1@ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)
