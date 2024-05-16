## Overview  
 * This repository contains scripts used for the analysis of single-cell RNAseq data from lung adenocarcinoma patients.
 * It contains `Rmarkdown` script in [Rmarkdown/current_pipeline](https://github.com/Yohan-Bosse-Lab/scRNA/tree/main/Rmarkdown/current_pipeline) than can be run in a `Rstudio` sesssion. Each contain specific analyses for this [preprint](https://www.biorxiv.org/content/10.1101/2024.02.20.581199v1.abstract).
     * `scRNA_qc_annotation.Rmd` (You need to run this first to QC the 10x data, annotate and prepare the Seurat objects)
     * `scRNA_NORMAL.Rmd` (Run this for Normal-specific analyses)
     * `scRNA_TUMOR.Rmd` (Run this for Tumor-specific analyses)
     * `scRNA_DEPLETED.Rmd` (Run this for Immune-depleted specific analyses)
     * `scRNA_ligand_receptor_interactome.Rmd` (cell-chat analyses)
     * `scRNA_annotation_plots.Rmd` (many plots are generated here)
     * `scRNA_umap_featureplots_allcells.Rmd` (two supplementary plots)
     * `scRNA_bulkDEG.Rmd` (differential gene expression analyses and plots)
     * `scRNA_H&E.Rmd` (immunohistochemical analyses and plots)
  * There are several helper functions described in the [R](https://github.com/Yohan-Bosse-Lab/scRNA/tree/main/R) subfolder directly.



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


## Further information
  * sebastien.renaut.1@ulaval.ca and yohan.bosse@criucpq.ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)


## Reference  
  *  Sébastien Renaut, Victoria Saavedra Armero, Dominique K. Boudreau, Nathalie Gaudreault, Patrice Desmeules, Sébastien Thériault, Patrick Mathieu, Philippe Joubert, Yohan Bossé. Single-cell and single-nucleus RNA-sequencing from paired normal-adenocarcinoma lung samples provides both common and discordant biological insights. *PLOS Genetics (accepted)*. 2024.
