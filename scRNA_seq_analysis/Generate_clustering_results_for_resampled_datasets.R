# Load libraries
library(Seurat)
library(tidyverse)
library(future)
library(harmony)
source("Hyperparameters_for_pipeline.R")
source("Pipelines.R")
plan("multicore", workers = number_of_thread)
# Record the start time
start_time <- Sys.time()
# load the reference atlas ccRCC data
scRNA_obj = readRDS(
  paste0(intermediate_dir,"columbia_british_individually_harmony_ccRCC_harmony.rds")
)

# bootstrap resampling the dataset with no replacement
rate = 0.8
ncells = ncol(scRNA_obj)
ncells.subsample <- round(ncells * rate)
scRNA_obj_resampled_list = list()
for(random_seed in seq(1,2,1)){
  set.seed(random_seed)
  subsample_index = sample(1:ncells, size = ncells.subsample, replace = F)
  scRNA_obj_resampled = CreateSeuratObject(counts = scRNA_obj[["RNA"]]@counts[,subsample_index])
  scRNA_obj_resampled = AddMetaData(scRNA_obj_resampled, 
                                    metadata = scRNA_obj@meta.data[rownames(scRNA_obj@meta.data)[subsample_index],])
  
  # do harmony batch correction pipeline on resampled datasets
  batch_key = "patient"
  scRNA_obj_resampled@meta.data[,batch_key] = as.character(factor(scRNA_obj_resampled@meta.data[,batch_key], 
                                                                  levels = unique(scRNA_obj_resampled@meta.data[,batch_key])))
  
  scRNA_obj_resampled = preprocess_standard_pipelines(scRNA_obj_resampled)
  scRNA_obj_resampled = RunPCA(scRNA_obj_resampled, npcs = n_of_pc)
  my_harmony_embeddings <- HarmonyMatrix(
    data_mat  = scRNA_obj_resampled@reductions$pca@cell.embeddings,
    meta_data = scRNA_obj_resampled@meta.data,
    vars_use  = c(batch_key),
    do_pca    = FALSE
  )
  #insert harmony reductions into Seurat object
  scRNA_obj_resampled@reductions$harmony = CreateDimReducObject(my_harmony_embeddings, assay = "RNA")
  #perform UMAP on harmony embeddings
  scRNA_obj_resampled <- RunUMAP(scRNA_obj_resampled, dims = 1:n_of_pc, reduction = "harmony", umap.method="umap-learn",
                                 metric="correlation", verbose = F, return.model = T)
  # find neighbours and clusters
  scRNA_obj_resampled = FindNeighbors(scRNA_obj_resampled, reduction = "harmony", dims = 1:n_of_pc)
  scRNA_obj_resampled = FindClusters(scRNA_obj_resampled, resolution = resolution_scope)
  scRNA_obj_resampled_list[[length(scRNA_obj_resampled_list)+1]] = scRNA_obj_resampled@meta.data
  end_time <- Sys.time()
  running_time = end_time - start_time
  write(paste0("Have done seed ", random_seed), stdout())
  write(paste0("The program has been running for ", running_time), stdout())
}

saveRDS(scRNA_obj_resampled_list, paste0(intermediate_dir, "columbia_british_individually_harmony_ccRCC_resampled_harmony_meta_df_100sample.rds"))