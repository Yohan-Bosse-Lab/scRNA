library(infercnv)


#a wrapper around infercnv functions to prepare data, run analysis and calculate a 'cnv score', something that, oddly, the infercnv library does not do.

infercnv_wrapper = function(seurat_obj = combined_permethod[[i]],
                            quantile = 0.5, 
                            out.dir = file.path(params$datapath,'infercnv/24samples/'),
                            level = 'predicted.ann_level_1',
                            gene_order_file = file.path(params$datapath,'gencode.v43.primary_assembly.annotation.positional')){
  #infercnv
  raw_counts_matrix = as.matrix(GetAssayData(object = seurat_obj, slot = "data"))
  quantile_threshold = quantile(rowSums(raw_counts_matrix),probs  = seq(0,1,by =0.1))[round(quantile*10)]
  raw_counts_matrix = raw_counts_matrix[rowSums(raw_counts_matrix) > quantile_threshold,]

  annotations_file = data.frame(annotations = seurat_obj@meta.data[,colnames(seurat_obj@meta.data) %in% level])
  rownames(annotations_file) = rownames(seurat_obj@meta.data)

  gene_order_file=read.delim(gene_order_file,header = F,row.names = 1)
  gene_order_file[,1] = as.factor(gene_order_file[,1])

  infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix, 
                                                          gene_order_file=gene_order_file,
                                                          annotations_file=annotations_file ,
                                                          ref_group_names=c("Immune"),
                                                          min_max_counts_per_cell = c(1, +Inf))

  infercnv_object_example <- infercnv::run(infercnv_object_example,
                                         cutoff=0.1,
                                         out_dir=out.dir,
                                         cluster_by_groups=T,
                                         denoise=TRUE,
                                         HMM=F,
                                         num_threads=12,
                                         no_plot=TRUE)


  ###cnv_Score (mean per gene)
  for(j in 1:nrow(annotations_file)){
    #cnv_mat[abs(infercnv_object_example@expr.data[,j]-1) > 0.05,j] = 1
    annotations_file$cnv_score[j] = mean((infercnv_object_example@expr.data[,j]-1)**2)
    }

  out_stats = annotations_file %>% 
    group_by(annotations) %>%
    summarise(n = mean(cnv_score),sd = sd(cnv_score),.groups='keep')

 # print(out_stats)

  seurat_obj@meta.data$cnv_score = annotations_file$cnv_score

  print(paste0('Done infercnvwrapper --- Time is: ',Sys.time())) 
  
  return(list(seurat_obj,out_stats))
}









