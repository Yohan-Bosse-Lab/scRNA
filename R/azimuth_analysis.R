#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)
library(arrow)


#this function is essentially just a wrapper around the Azimuth Shiny App.
azimuth_annotation = function(query = NULL,annotation_threshold = 0,
                              reference_path = "C:/Users/renseb01/Documents/scRNA/scRNA/data/HLCA/lung_2.0.0/",
                              v2=T,
                              version2.0.1="C:/Users/renseb01/Documents/scRNA/scRNA/data/HLCA/lung_2.0.1/data/annotations.parquet") {

  if(class(query)!='Seurat') stop("You did not provide an object of class Seurat")
  
  input_data = query
  
# Load the reference
reference <- LoadReference(path = reference_path)

#new reference v2.0.1
if(v2) {
  annotations <- read_parquet(version2.0.1)
  annotations = annotations[annotations$`__index_level_0__` %in% rownames(reference$map@meta.data),]
  
  reference$map@meta.data$ann_level_1 = annotations$ann_level_1
  reference$map@meta.data$ann_level_2 = annotations$ann_level_2
  reference$map@meta.data$ann_level_3 = annotations$ann_level_3
  reference$map@meta.data$ann_level_4 = annotations$ann_level_4
  reference$map@meta.data$ann_level_5 = annotations$ann_level_5
  reference$map@meta.data$ann_finest_level = annotations$ann_finest_level
  
  print('re-annotating reference with version 2.0.1')
}

# Load the query object for mapping
# Change the file path based on where the query file is located on your system.
#query <- LoadFileInput(path = "C:/Users/renseb01/Documents/scRNA/scRNA/data/filtered_feature_bc_matrix.h5")
query <- ConvertGeneNames(
  object = query,
  reference.names = rownames(x = reference$map),
  homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
)

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_level_5","ann_finest_level"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_level_5", "ann_finest_level")
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

NNTransform <- function(
    object,
    meta.data,
    neighbor.slot = "query_ref.nn",
    key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)
#get the cluster preservation score:
qps = ClusterPreservationScore(query,5000)
print(qps)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

new_metadata_columns = query@meta.data[,regexpr('ann_|mapping.score',colnames(query@meta.data))>0]

for(i in 1:6)
  {
  new_metadata_columns[new_metadata_columns[,(i*2)-1]<annotation_threshold,(i*2)] = 'unknown'
  new_metadata_columns[new_metadata_columns[,(i*2)-1]<annotation_threshold,(i*2)] = 'unknown'
  
  print(paste0(colnames(new_metadata_columns)[(i*2)],' --- ',
               length(new_metadata_columns[new_metadata_columns[,(i*2)-1]<annotation_threshold,(i*2)]),' values out of ',nrow(new_metadata_columns),' set to unknown'))
  }
  


###
predicted_levels = c('predicted.ann_finest_level','predicted.ann_level_1','predicted.ann_level_2','predicted.ann_level_3','predicted.ann_level_4','predicted.ann_level_5')
scores_levels = c('predicted.ann_finest_level.score','predicted.ann_level_1.score','predicted.ann_level_2.score','predicted.ann_level_3.score','predicted.ann_level_4.score','predicted.ann_level_5.score')


for(i in 1:6)
{
  predicted = query@meta.data[,colnames(query@meta.data) %in% predicted_levels[i]]
  scores = query@meta.data[,colnames(query@meta.data) %in% scores_levels[i]]

  predicted[scores<annotation_threshold] = 'unknown'
  
  query@meta.data[,colnames(query@meta.data) %in% predicted_levels[i]] = predicted
  
}

new_metadata_columns

input_data@meta.data = cbind(input_data@meta.data,new_metadata_columns)

return(list(input_data,qps,query))

}


#cluster preservation score from Azimuth 
ClusterPreservationScore <- function(query, ds.amount) {
  query <- DietSeurat(object = query, assays = "refAssay", scale.data = TRUE, counts = FALSE, dimreducs = "integrated_dr")
  if (ncol(x = query) > ds.amount) {
    query <- subset(x = query, cells = sample(x = Cells(x = query), size = ds.amount))
  }
  dims <- min(50, getOption(x = "Azimuth.map.ndims"))
  query <- RunPCA(object = query, npcs = dims, verbose = FALSE)
  query <- FindNeighbors(
    object = query,
    reduction = 'pca',
    dims = 1:dims,
    graph.name = paste0("pca_", c("nn", "snn"))
  )
  query[["orig_neighbors"]] <- as.Neighbor(x = query[["pca_nn"]])
  query <- FindClusters(object = query, resolution = 0.6, graph.name = 'pca_snn')
  query <- FindNeighbors(
    object = query,
    reduction = 'integrated_dr',
    dims = 1:dims,
    return.neighbor = TRUE,
    graph.name ="integrated_neighbors_nn"
  )
  ids <- Idents(object = query)
  integrated.neighbor.indices <- Indices(object = query[["integrated_neighbors_nn"]])
  proj_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- integrated.neighbor.indices[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = proj_ent) <- Cells(x = query)
  orig.neighbor.indices <- Indices(object = query[["orig_neighbors"]])
  orig_ent <- unlist(x = lapply(X = 1:length(x = Cells(x = query)), function(x) {
    neighbors <- orig.neighbor.indices[x, ]
    nn_ids <- ids[neighbors]
    p_x <- prop.table(x = table(nn_ids))
    nn_entropy <- sum(p_x * log(x = p_x), na.rm = TRUE)
    return(nn_entropy)
  }))
  names(x = orig_ent) <- Cells(x = query)
  stat <- median(
    x = tapply(X = orig_ent, INDEX = ids, FUN = mean) -
      tapply(X = proj_ent, INDEX = ids, FUN = mean)
  )
  if (stat <= 0) {
    stat <- 5.00
  } else {
    stat <- -1 * log2(x = stat)
    stat <- MinMax(data = stat, min = 0.00, max = 5.00)
  }
  return(stat)
}



#relabel all unknowns as unclassified in a data.frame or seurat object.
relabeler = function(seurat.data=data){

  levels = c('predicted.ann_level_1','predicted.ann_level_2','predicted.ann_level_3','predicted.ann_level_5','predicted.ann_level_5','predicted.ann_finest_level')
  
  if(class(seurat.data) == 'Seurat') {
    seurat.data@meta.data$predicted.ann_level_1[seurat.data@meta.data$predicted.ann_level_1 =='unknown'] = 'Unclassified'
    seurat.data@meta.data$predicted.ann_level_2[seurat.data@meta.data$predicted.ann_level_2 =='unknown'] = 'Unclassified'
    seurat.data@meta.data$predicted.ann_level_3[seurat.data@meta.data$predicted.ann_level_3 =='unknown'] = 'Unclassified'
    seurat.data@meta.data$predicted.ann_level_4[seurat.data@meta.data$predicted.ann_level_4 =='unknown'] = 'Unclassified'
    seurat.data@meta.data$predicted.ann_level_5[seurat.data@meta.data$predicted.ann_level_5 =='unknown'] = 'Unclassified'
    seurat.data@meta.data$predicted.ann_finest_level[seurat.data@meta.data$predicted.ann_finest_level=='unknown'] = 'Unclassified'
  }
  
  
  if(class(seurat.data) == 'data.frame') {
    for(i in c(1:6)[levels %in% colnames(seurat.data)]){
      seurat.data[,colnames(seurat.data) == levels[i]][seurat.data[,colnames(seurat.data) == levels[i]] == 'unknown'] = 'Unclassified'
    }
  }
  
  if(class(seurat.data) == 'list'){
    for(i in 1:length(seurat.data)){
    
      if(length(seurat.data[[i]])>1){
        for(j in 1:length(seurat.data[[i]])){
          seurat.data[[i]][[j]]@meta.data$predicted.ann_level_1[seurat.data[[i]][[j]]@meta.data$predicted.ann_level_1=='unknown'] = 'Unclassified'
          seurat.data[[i]][[j]]@meta.data$predicted.ann_level_2[seurat.data[[i]][[j]]@meta.data$predicted.ann_level_2=='unknown'] = 'Unclassified'
          seurat.data[[i]][[j]]@meta.data$predicted.ann_level_3[seurat.data[[i]][[j]]@meta.data$predicted.ann_level_3=='unknown'] = 'Unclassified'
          seurat.data[[i]][[j]]@meta.data$predicted.ann_level_4[seurat.data[[i]][[j]]@meta.data$predicted.ann_level_4=='unknown'] = 'Unclassified'
          seurat.data[[i]][[j]]@meta.data$predicted.ann_level_5[seurat.data[[i]][[j]]@meta.data$predicted.ann_level_5=='unknown'] = 'Unclassified'
          seurat.data[[i]][[j]]@meta.data$predicted.ann_finest_level[seurat.data[[i]][[j]]@meta.data$predicted.ann_finest_level=='unknown'] = 'Unclassified'
        }
      }
    
      if(length(seurat.data[[i]])==1){
        seurat.data[[i]]@meta.data$predicted.ann_level_1[seurat.data[[i]]@meta.data$predicted.ann_level_1=='unknown'] = 'Unclassified'
        seurat.data[[i]]@meta.data$predicted.ann_level_2[seurat.data[[i]]@meta.data$predicted.ann_level_2=='unknown'] = 'Unclassified'
        seurat.data[[i]]@meta.data$predicted.ann_level_3[seurat.data[[i]]@meta.data$predicted.ann_level_3=='unknown'] = 'Unclassified'
        seurat.data[[i]]@meta.data$predicted.ann_level_4[seurat.data[[i]]@meta.data$predicted.ann_level_4=='unknown'] = 'Unclassified'
        seurat.data[[i]]@meta.data$predicted.ann_level_5[seurat.data[[i]]@meta.data$predicted.ann_level_5=='unknown'] = 'Unclassified'
        seurat.data[[i]]@meta.data$predicted.ann_finest_level[seurat.data[[i]]@meta.data$predicted.ann_finest_level=='unknown'] = 'Unclassified'
      }
    }
  }
  return(seurat.data) 
}



