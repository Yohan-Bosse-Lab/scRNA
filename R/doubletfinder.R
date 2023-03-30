library(DoubletFinder)
library(Seurat)


# A wrapper for DoubletFinder
doubletfinder_wrapper = function(seurat.object= NULL,pK_empirical=NULL,db_rate = 0.05,pN = 0.25,PCs = 1:30){
  a=Sys.time()
  if(is.null(pK_empirical)){
    sweep.res.list <- paramSweep_v3(seurat.object, PCs = PCs, sct = T)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_empirical = as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]))
  }
  
  ## Run DoubletFinder 
  out <- doubletFinder_v3(seurat.object,
                               PCs = PCs,
                               pN = pN,
                               pK = pK_empirical,
                               nExp = round(db_rate*nrow(seurat.object@meta.data)),
                               reuse.pANN = FALSE,
                               sct = TRUE)
  
  #rename column in order to use it below to subset
  colnames(out@meta.data)[regexpr('DF.classification',colnames(out@meta.data),fixed =T)>0]= 'DF.classifications'
  
  # Subset dataset
  out <- subset(out, subset = DF.classifications == 'Singlet')

  print(paste0('Empirical pk is : ',pK_empirical))
  print(Sys.time()-a)
  return(list(out,pK_empirical))
}
