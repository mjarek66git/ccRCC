set.seed(10)
library(Seurat)
library(future)
library(harmony)
source("Azimuth.R")
source("Pipelines.R")
source("Hyperparameters_for_pipeline.R")
source("tools.R")
check_dir(intermediate_dir)
options(future.globals.maxSize= 80*1000*1024^2)
plan("multicore", workers = number_of_thread)
# load data from the file "raw_data.rds"
# "raw_data.rds" can be downloaded using this link:https://drive.google.com/file/d/1CEdUnspMfthex0XSxemxCy1j5TWkWxCE/view?usp=sharing
scRNA_objs = readRDS("your_path_to_raw_data.rds")
for(i in seq(1,2,1)){
  scRNA_obj = scRNA_objs[[i]]
  study_source = unique(scRNA_obj@meta.data[,"study"])
  min_cells_expressed_per_gene = 20
  if(study_source == "columbia"){
    selected_f <- Matrix::rowSums(x = scRNA_obj[["RNA"]]@counts > 0) >= min_cells_expressed_per_gene
    scRNA_obj = scRNA_obj[selected_f,]
    scRNA_obj = PercentageFeatureSet(scRNA_obj, pattern = "^MT-", col.name = "percent.mt")
    png("QC.png", width = 16, height = 9, res = 300, units = "in")
    VlnPlot(scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
    dev.off()
    nCount_RNA_min = 1000
    nCount_RNA_max = 30000
    nFeature_RNA_min = 500
    nFeature_RNA_max = 5000
    patient_key = "patient"
  }else if(study_source == "british"){
    Idents(scRNA_obj) = "british"
    selected_f <- Matrix::rowSums(x = scRNA_obj[["RNA"]]@counts > 0) >= min_cells_expressed_per_gene
    scRNA_obj = scRNA_obj[selected_f,]
    scRNA_obj = PercentageFeatureSet(scRNA_obj, pattern = "^MT-", col.name = "percent.mt")
    png("QC.png", width = 16, height = 9, res = 300, units = "in")
    VlnPlot(scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
    dev.off()
    nCount_RNA_min = 1000
    nCount_RNA_max = 6000
    nFeature_RNA_min = 500
    nFeature_RNA_max = 3000
    patient_key = "patient"
  }
  scRNA_obj = subset(scRNA_obj, subset = percent.mt <= mt_threshold & 
                       nCount_RNA >= nCount_RNA_min & nCount_RNA <= nCount_RNA_max &
                       nFeature_RNA >= nFeature_RNA_min & nFeature_RNA <= nFeature_RNA_max)
  batch_correction_method = batch_correction_methods[3]
  scRNA_obj = preprocess_standard_pipelines(scRNA_obj)
  scRNA_obj = RunPCA(scRNA_obj, npcs = n_of_pc)
  my_harmony_embeddings <- HarmonyMatrix(
    data_mat  = scRNA_obj@reductions$pca@cell.embeddings,
    meta_data = scRNA_obj@meta.data,
    vars_use  = patient_key,
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
  scRNA_obj_azimuth = run_azimuth(scRNA_obj, reference_path = "./azimuth_reference/")
  scRNA_obj@meta.data[,c("azimuth_anno_l1_score", "azimuth_anno_l1", "azimuth_anno_l2_score", 
                         "azimuth_anno_l2", "azimuth_anno_l3_score", "azimuth_anno_l3")]=
    scRNA_obj_azimuth@meta.data[,(ncol(scRNA_obj_azimuth@meta.data)-5):ncol(scRNA_obj_azimuth@meta.data)]
  rm(scRNA_obj_azimuth)
  saveRDS(scRNA_obj, paste0(intermediate_dir, paste(study_source, batch_correction_method, n_of_pc, 
                                                    "PCs", "azimuthed.rds", sep="_")))
}