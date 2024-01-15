preprocess_standard_pipelines = function(obj, n_variable_features = 3000){
  # Normalize data
  obj = NormalizeData(obj, verbose = F)
  # FindVariableFeatures
  obj = FindVariableFeatures(obj, nfeature = n_variable_features, verbose = F)
  # Scale data
  obj = ScaleData(obj, verbose = F)
  
  return(obj)
}

preprocess_SCT = function(obj, n_variable_features = 3000){
  obj = SCTransform(obj, variable.features.n = n_variable_features, 
                    vars.to.regress = "percent.mt", verbose = F)
  return(obj)
}

run_pca_and_umap = function(obj, n_of_pc = 30){
  # Run PCA
  obj <- RunPCA(obj, npcs = n_of_pc, verbose = F)
  # Run Umap
  obj <- RunUMAP(obj, dims = 1:n_of_pc, reduction = "pca", umap.method="umap-learn",
                 metric="correlation", verbose = F, return.model = T)
  return(obj)
}

find_neighbours_and_clusters = function(obj, n_of_pc = 30, res = seq(0.1,1,0.1)){
  obj <- FindNeighbors(obj, dims = 1:n_of_pc, verbose = F)
  obj <- FindClusters(obj, resolution = res)
  return(obj)
}

options(future.globals.maxSize = 20000 * 1024^2)
#function to integrate a splited seurat datasets with split key and 
#ref as the dataset who has the largest number of cells
Integrate_seurat_obj_by_key = function(obj, key, method_SCT = T, method_rpca = F,
                                       n.features=2000, n.anchors=5, k.weight = 100){
  obj_list = SplitObject(obj, split.by = key)
  if(method_SCT){
    if(method_rpca){
      obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi", variable.features.n = n.features)
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
      obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", 
                                        anchor.features = features, reduction = "rpca",
                                        k.anchor = n.anchors)
    }else{
      obj_list <- lapply(X = obj_list, FUN = SCTransform, variable.features.n = n.features)
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", 
                                        anchor.features = features,
                                        k.anchor = n.anchors)
    }
    result_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F,
                                 k.weight = k.weight)
  }else{
    
    if(method_rpca){
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n.features)
      })
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, 
                                        reduction = "rpca", 
                                        k.anchor = n.anchors)
    }else{
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n.features)
      })
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, 
                                        k.anchor = n.anchors)
    }
    result_data <- IntegrateData(anchorset = anchors, verbose = F, k.weight = k.weight)
    DefaultAssay(result_data) <- "integrated"
    result_data <- ScaleData(result_data, verbose = FALSE)
  }
  
  return(result_data)
}

