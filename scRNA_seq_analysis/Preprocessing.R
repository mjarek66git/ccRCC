set.seed(10)
library(Seurat)
CreateSeuratObject()
library(future)
library(harmony)
source("Azimuth.R")
source("Pipelines.R")
source("Hyperparameters_for_pipeline.R")
options(future.globals.maxSize= 80*1000*1024^2)
plan("multicore", workers = number_of_thread)
#load data from all cohort list
scRNA_objs = readRDS("/scratch/yanv5j/all_cohorts_in_list.rds")
for(i in seq(1,4,1)){
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
  }else if(study_source == "harvard"){
    Idents(scRNA_obj) = "harvard"
    selected_f <- Matrix::rowSums(x = scRNA_obj[["RNA"]]@counts > 0) >= min_cells_expressed_per_gene
    scRNA_obj = scRNA_obj[selected_f,]
    scRNA_obj = PercentageFeatureSet(scRNA_obj, pattern = "^MT-", col.name = "percent.mt")
    png("QC.png", width = 16, height = 9, res = 300, units = "in")
    VlnPlot(scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
    dev.off()
    nCount_RNA_min = 1000
    nCount_RNA_max = 25000
    nFeature_RNA_min = 500
    nFeature_RNA_max = 5000
    patient_key = "patient"
  }else if(study_source == "michigan"){
    Idents(scRNA_obj) = "michigan"
    selected_f <- Matrix::rowSums(x = scRNA_obj[["RNA"]]@counts > 0) >= min_cells_expressed_per_gene
    scRNA_obj = scRNA_obj[selected_f,]
    scRNA_obj = PercentageFeatureSet(scRNA_obj, pattern = "^MT-", col.name = "percent.mt")
    png("QC.png", width = 16, height = 9, res = 300, units = "in")
    VlnPlot(scRNA_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
    dev.off()
    nCount_RNA_min = 1000
    nCount_RNA_max = 25000
    nFeature_RNA_min = 500
    nFeature_RNA_max = 6000
    patient_key = "patient"
  }
  scRNA_obj = subset(scRNA_obj, subset = percent.mt <= mt_threshold & 
                       nCount_RNA >= nCount_RNA_min & nCount_RNA <= nCount_RNA_max &
                       nFeature_RNA >= nFeature_RNA_min & nFeature_RNA <= nFeature_RNA_max)
  batch_correction_method = batch_correction_methods[3]
  if(batch_correction_method == "harmony"){
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
  }else if(batch_correction_method == "Seurat"){
    scRNA_obj = Integrate_seurat_obj_by_key(
      scRNA_obj, key = patient_key,
      method_SCT = F, method_rpca = F
    )
    scRNA_obj = run_pca_and_umap(scRNA_obj, n_of_pc = n_of_pc)
    scRNA_obj = find_neighbours_and_clusters(scRNA_obj, n_of_pc = n_of_pc, res = resolution_scope)
  }else if(batch_correction_method == "SeuratSCT"){
    scRNA_obj = Integrate_seurat_obj_by_key(
      scRNA_obj, key = patient_key,
      method_SCT = T, method_rpca = F
    )
    scRNA_obj = run_pca_and_umap(scRNA_obj, n_of_pc = n_of_pc)
    scRNA_obj = find_neighbours_and_clusters(scRNA_obj, n_of_pc = n_of_pc, res = resolution_scope)
  }else if(batch_correction_method == "scVI"){
    #import additional libraries
    library(reticulate)
    sc <- import('scanpy', convert = FALSE)
    scvi <- import('scvi', convert = FALSE)
    scvi$settings$progress_bar_style = 'tqdm'
    
    # normalize data
    scRNA_obj = NormalizeData(scRNA_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    # find highly variable genes
    scRNA_obj = FindVariableFeatures(scRNA_obj, selection.method = "vst", nfeatures = n_of_highly_variable_genes)
    top2000 = head(VariableFeatures(scRNA_obj), n_of_highly_variable_genes)
    # initialize ann data
    scRNA_adata <- sc$AnnData(
      X   = t(as.matrix(GetAssayData(scRNA_obj[top2000], slot='counts'))),
      obs = scRNA_obj[top2000][[]]
    )
    # run setup_anndata, use column patient for batch
    scvi$model$SCVI$setup_anndata(scRNA_adata, batch_key = 'patient')
    # create the model
    model = scvi$model$SCVI(scRNA_adata, n_latent=as.integer(n_of_pc), n_layers=as.integer(2))
    # transfer model to GPU devices
    model$to_device('cuda:0')
    # train the model
    model$train(max_epochs = as.integer(250))
    # get latent representation
    latent = model$get_latent_representation()
    # put it back in our original Seurat object
    latent <- as.matrix(latent)
    rownames(latent) = colnames(scRNA_obj)
    scRNA_obj[['scvi']] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(scRNA_obj))
    scRNA_obj <- RunUMAP(scRNA_obj, dims = 1:n_of_pc, reduction = "scvi", umap.method="umap-learn",
                         metric="correlation", verbose = F, return.model = T)
    scRNA_obj <- FindNeighbors(scRNA_obj, dims = 1:n_of_pc, reduction = "scvi")
    scRNA_obj <- FindClusters(scRNA_obj, resolution = resolution_scope)
  }
  scRNA_obj_azimuth = run_azimuth(scRNA_obj, reference_path = "./azimuth_reference/")
  scRNA_obj@meta.data[,c("azimuth_anno_l1_score", "azimuth_anno_l1", "azimuth_anno_l2_score", 
                         "azimuth_anno_l2", "azimuth_anno_l3_score", "azimuth_anno_l3")]=
    scRNA_obj_azimuth@meta.data[,(ncol(scRNA_obj_azimuth@meta.data)-5):ncol(scRNA_obj_azimuth@meta.data)]
  rm(scRNA_obj_azimuth)
  saveRDS(scRNA_obj, paste0("../scRNA_data/", paste(study_source, batch_correction_method, n_of_pc, 
                                                    "PCs", "azimuthed.rds", sep="_")))
}