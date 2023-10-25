## Overview  
 * This repository contains scripts used for the analysis of single-cell RNAseq data from lung adenocarcinoma patients.
 * It contains`Rmarkdown` script than can be run as such from the command line (or alternatively, `Knit` directly in `R/Rstudio`) :
 * It describes the analyses I did for this manuscript (in preparation): *Single-cell and single-nucleus RNA-sequencing of the same paired normal-adenocarcinoma lung samples provides divergent biological information*
     * `scRNA_qc_annotation.Rmd` (You need to run this first to QC the 10x data, annotate and prepare the Seurat objects)
     * `scRNA_NORMAL.Rmd` (Run this for Normal-specific analyses)
     * `scRNA_TUMOR.Rmd` (Run this for Tumor-specific analyses)
     * `scRNA_DEPLETED.Rmd` (Run this for Immune-depleted specific analyses)
     * `ligand_receptor_interactome.Rmd` (cell-chat analyses)
     * `scRNA_annotation_plots.Rmd` (many plots are generated here)
     * `scRNA_umap_featureplots_allcells.Rmd` (two supplementary plots)


## Installation
  * This is R-based, libraries requires are described in `.Rmd` files, but mainly:
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
Rscript -e "rmarkdown::render('scRNA_qc_annotation.Rmd',params = list(nFeature_min=250, 
nCount_min=500, nCount_max=50000, percent_MT_max=3)"
```

## Further information
  * sebastien.renaut.1@ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)
