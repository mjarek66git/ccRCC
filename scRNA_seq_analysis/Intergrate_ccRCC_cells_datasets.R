library(Seurat)
source("Pipelines.R")
source("Hyperparameters_for_pipeline.R")
source("tools.R")
check_dir(intermediate_dir)
ATP_no_F = readRDS("ATP_dic.rds")
study_indexes = seq(1,2,1)
scRNA_obj_list = list()
for(i in study_indexes){
  scRNA_obj_individual = readRDS(paste0(intermediate_dir, study_sources[i], "_harmony_ccRCC.rds"))
  # correction for ATP genes
  if(any(rownames(scRNA_obj_individual[["RNA"]]@counts) %in% ATP_no_F)){
    rownames(scRNA_obj_individual[["RNA"]]@counts)[rownames(scRNA_obj_individual[["RNA"]]@counts) %in% ATP_no_F] = names(ATP_no_F)
    rownames(scRNA_obj_individual[["RNA"]]@data)[rownames(scRNA_obj_individual[["RNA"]]@data) %in% ATP_no_F] = names(ATP_no_F)
  }
  scRNA_obj_list[[i]] = scRNA_obj_individual
}
scRNA_obj = merge(scRNA_obj_list[[1]], scRNA_obj_list[2:length(scRNA_obj_list)])
# remove some non-expressed features
selected_f <- Matrix::rowSums(x = scRNA_obj[["RNA"]]@counts > 0) >= 10
scRNA_obj = scRNA_obj[selected_f,]
# define the batch key
batch_key = "patient"
# run integration
library(future)
plan("multicore", workers = number_of_thread)
batch_correction_method = batch_correction_methods[3]
print(paste0("batch correction method is ", batch_correction_method))
scRNA_obj@meta.data[,batch_key] = as.character(factor(scRNA_obj@meta.data[,batch_key], 
                                                      levels = unique(scRNA_obj@meta.data[,batch_key])))
if(batch_correction_method == "harmony"){
  library(harmony)
  scRNA_obj = preprocess_standard_pipelines(scRNA_obj)
  scRNA_obj = RunPCA(scRNA_obj, npcs = n_of_pc)
  my_harmony_embeddings <- HarmonyMatrix(
    data_mat  = scRNA_obj@reductions$pca@cell.embeddings,
    meta_data = scRNA_obj@meta.data,
    vars_use  = batch_key,
    do_pca    = FALSE
  )
  #insert harmony reductions into Seurat object
  scRNA_obj@reductions$harmony = CreateDimReducObject(my_harmony_embeddings, assay = "RNA")
  #perform UMAP on harmony embeddings
  scRNA_obj <- RunUMAP(scRNA_obj, dims = 1:n_of_pc, reduction = "harmony", umap.method="umap-learn",
                       metric="correlation", verbose = F, return.model = T)
  # find neighbours and clusters
  scRNA_obj = FindNeighbors(scRNA_obj, reduction = "harmony", dims = 1:n_of_pc)
  scRNA_obj = FindClusters(scRNA_obj, resolution = resolution_scope)
}

saveRDS(scRNA_obj, paste0(intermediate_dir, "columbia_british_individually_harmony_ccRCC_harmony.rds"))