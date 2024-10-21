import anndata
import pandas as pd

#(This is done on GENNO, otherwise I don't have enough CPU on laptop)

#load data
bigdata = anndata.read_h5ad('hlca_dataset_2.8M_cells.h5ad',backed= 'r')

#diseases of interest
disease_of_interest = ["interstitial lung disease","lung adenocarcinoma","lung large cell carcinoma","pleomorphic carcinoma","pulmonary sarcoidosis","lung large cell carcinoma","squamous cell lung carcinoma"]
#disease_of_interest = ["lung large cell carcinoma"]

#subset
bigdata_subset = bigdata[bigdata.obs.disease.isin(disease_of_interest)]

#write
bigdata_subset.write('hlca_subset_2.8M_cells.h5ad',compression = 'gzip')

metadata = bigdata_subset.obs 
metadata.to_csv('hlca_subset_2.8M_cells.csv')


###
###
###
#This is done in R (on GENNO)
fix_convert = function(h5seurat = "bigdata_subset.h5seurat") {

  f <- hdf5r::H5File$new(h5seurat, "r+")
  groups <- f$ls(recursive = TRUE)
  
  #
  print(paste0('The number of categories is: ',length(groups$name[grepl("categories", groups$name)])))
  for (name in groups$name[grepl("categories", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "levels")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
  }
  
  #
  print(paste0('The number of codes is: ',length(groups$name[grepl("codes", groups$name)])))
  for (name in groups$name[grepl("codes", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "values")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
    grp <- f[[new_name]]
    grp$write(args = list(1:grp$dims), value = grp$read() + 1)
  }

    f$close_all()

    return('fixed the levels & groups conversion')
  }

#
SeuratDisk::Convert("hlca_subset_2.8M_cells.h5ad", dest = "h5seurat", overwrite = TRUE)
fix_convert('hlca_subset_2.8M_cells.h5seurat')


#####
# Loading is done in R (on laptop, otherwise it errors out with old version of SeuratDisk...)
#####
bigdata_subset <- SeuratDisk::LoadH5Seurat("hlca_subset_2.8M_cells.h5seurat",meta.data = FALSE,misc=FALSE)
metadata = read.csv('hlca_subset_2.8M_cells.csv',row.names =  1)

bigdata_subset@meta.data = metadata