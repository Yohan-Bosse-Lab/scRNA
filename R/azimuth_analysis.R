#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)
library(arrow)


azimuth_annotation = function(query = NULL,annotation_threshold = 0, reference_path = "C:/Users/renseb01/Documents/scRNA/scRNA/data/HLCA/lung_2.0.0/",v2=T,version2.0.1="C:/Users/renseb01/Documents/scRNA/scRNA/data/HLCA/lung_2.0.1/data/annotations.parquet") {

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

  predicted[scores<0.5] = 'unknown'
  
  query@meta.data[,colnames(query@meta.data) %in% predicted_levels[i]] = predicted
  
}



new_metadata_columns


input_data@meta.data = cbind(input_data@meta.data,new_metadata_columns)

return(list(input_data,qps,query))

}
# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
#id <- c("ann_level_1","ann_level_2","ann_level_3","ann_level_4", "ann_finest_level")
#predicted.id <- paste0("predicted.", id)

# DimPlot of the reference
#DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
#DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id[2], label = TRUE) + NoLegend()





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



if(1==2){
#POC

#setup
packageVersion(pkg = "Seurat")
#[1] ‘4.3.0’
packageVersion(pkg = "Azimuth")
#[1] ‘0.4.6’

SeuratData::InstallData('pbmc3k')
pbmc3k = SeuratData::LoadData('pbmc3k')
SeuratDisk::SaveH5Seurat(pbmc3k,'pbmc3k.h5Seurat')

#run Azimuth locally 
pbmc3k = RunAzimuth(pbmc3k,reference = 'pbmcref')

#run through Shiny App
azimuth_pred_App = predictions <- read.delim('azimuth_pred.tsv', row.names = 1)

#make sure we have the same cells
all.equal(rownames(azimuth_pred_App),rownames(pbmc3k@meta.data))

#plot
plot(pbmc3k@meta.data$mapping.score,azimuth_pred_App$mapping.score)

}



if(1==2){
#taken from the https://github.com/satijalab/azimuth-references/tree/master/human_lung_v2
#prep reference
library(Seurat)
library(SeuratObject)
library(Azimuth)
library(Matrix)
library(arrow)

setwd('scRNA/scRNA/data/HLCA/lung_2.0.1/')


#script = "scripts/build_reference.R",
counts.path  = "data/counts.rds"
annotations.path = "data/annotations.parquet"
dr.path = "data/scanvi.parquet"
ref.path = "reference/ref.Rds"
annoy.path  = "reference/idx.annoy"
full.obj.path  ="full_reference.Rds"

mtx <- readRDS(counts.path)
obj <- CreateSeuratObject(counts = mtx)

# load annotations
annotations <- read_parquet(annotations.path)
annotations <- as.data.frame(annotations)[,c("ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level")]
rownames(annotations) <- Cells(obj)
obj <- AddMetaData(obj, metadata = annotations)
print(head(obj))

# load in the scANVI latent space which contains 30 dimensions
latent.space <- as.matrix(read_parquet(dr.path))
rownames(latent.space) <- Cells(obj)
scanvi.dr <- CreateDimReducObject(embeddings = as.matrix(latent.space), key = "SCANVI")
obj[["scanvi"]] <- scanvi.dr

# find neighbors based on scANVI latent space
obj <- FindNeighbors(obj, reduction = "scanvi")

# Run SCTransform on the raw counts

###DOES NOT WORK
#HERE HERE HERE 
#NOT ENOUGH MEMORE V2 BREKS ...
obj <- SCTransform(obj,method = glmGamPoi_offset, n_cells=2000,exclude_poisson = TRUE, vst.flavor = "v2")

# run sPCA
obj <- RunSPCA(object = obj, assay = "SCT", graph = "RNA_snn")

# Force RunUMAP to run with n.epochs to prevent RunUMAP from running for 0 epochs
# Related: https://scanpy.discourse.group/t/umap-incorrectly-installed/663
obj <- RunUMAP(obj, dims = 1:30, reduction = "scanvi", n.epochs = 200, return.model = TRUE)

# save the full size object to perform marker identification
saveRDS(obj, file = full.obj.path)

annotations <- c("ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level")
ref <- AzimuthReference(obj,
                        refUMAP = 'umap',
                        refDR = 'spca',
                        dims = 1:50, # use 50 dimensions from the sPCA dimensional reduction
                        plotref = 'umap',
                        reference.version = '2.0.1',
                        metadata = annotations)

saveRDS(object = ref, file = ref.path, compress = F)
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]],
               file = annoy.path)


}