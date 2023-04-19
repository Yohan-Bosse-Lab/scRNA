#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)


azimuth_annotation = function(query = NULL,reference_path = "C:/Users/renseb01/Documents/scRNA/scRNA/data/HLCA/data_azimuth/") {

  if(class(query)!='Seurat') stop("You did not provide an object of class Seurat")
  
  input_data = query
  
# Load the reference
reference <- LoadReference(path = reference_path)

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

refdata <- lapply(X = c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_finest_level"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- c("ann_level_1","ann_level_2","ann_level_3","ann_level_4", "ann_finest_level")
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

for(i in 1:5)
  {
  new_metadata_columns[new_metadata_columns[,(i*2)-1]<0.4,(i*2)] = 'unknown'
  
  print(paste0(colnames(new_metadata_columns)[(i*2)],' --- ',
               length(new_metadata_columns[new_metadata_columns[,(i*2)-1]<0.4,(i*2)]),' values out of ',nrow(new_metadata_columns),' set to unknown'))
  }
  

  

new_metadata_columns


input_data@meta.data = cbind(input_data@meta.data,new_metadata_columns)

return(input_data)

}
# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
#id <- c("ann_level_1","ann_level_2","ann_level_3","ann_level_4", "ann_finest_level")
#predicted.id <- paste0("predicted.", id)

# DimPlot of the reference
#DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
#DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id[2], label = TRUE) + NoLegend()

