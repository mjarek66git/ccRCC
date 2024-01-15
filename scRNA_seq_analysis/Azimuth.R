library(Azimuth)

run_azimuth <- function(seurat_obj, reference_path = "/Volumes/Jason/docker_dir/references/"){
  
  #load the reference data
  reference <- LoadReference(path = reference_path)
  
  #save rds
  saveRDS(seurat_obj, "azimuth.rds")
  
  #load query
  query <- LoadFileInput(path = "azimuth.rds")
  
  #Convert gene names
  query <- ConvertGeneNames(
    object = query,
    reference.names = rownames(x = reference$map),
    homolog.table = './homologs.rds'
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
    dims = 1:100,
    n.trees = 20,
    mapping.score.k = 100
  )
  
  # Transfer cell type labels and impute protein expression
  refdata <- lapply(X = c("annotation.l1", "annotation.l2", "annotation.l3"), function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- c("annotation.l1", "annotation.l2", "annotation.l3")
  
  query <- TransferData(
    reference = reference$map,
    query = query,
    dims = 1:100,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE
  )
  
  # Calculate the embeddings of the query data on the reference
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE
  )
  
  return(query)
}
