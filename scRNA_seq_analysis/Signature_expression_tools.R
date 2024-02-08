#define colnames for different algorithms of calculating super genes expression
module_score_tag = "_ModuleScore"
z_score_tag = "_ZScore"
normalized_count_tag = "_NormalizedCount"
raw_count_tag = "_RawCount"
# define function to add z-score count of super gene to object
AddScaledNormalizeCount = function(object, gene_list){
  mat = object[["RNA"]]@data
  result_list = list()
  for(x in 1:length(gene_list)){
    selected_genes = gene_list[[x]]
    selected_genes = trimws(selected_genes, which = c("both"))
    marker_tag = names(gene_list[x])
    is_in_mat = selected_genes %in% rownames(mat)
    if(!all(is_in_mat)){
      not_available_genes = selected_genes[!is_in_mat]
      print(paste0(paste(not_available_genes, collapse = ","), 
                   " is not available for analysis"))
    }
    selected_genes = selected_genes[is_in_mat]
    selected_mat = mat[selected_genes, ]
    composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
    composite_marker_exp = scale(composite_marker_exp)
    marker_tag = paste0(marker_tag, z_score_tag)
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}
AddScaledNormalizeCountPerStage = function(object, gene_list){
  result_list = c()
  for(x in c("pT1b","pT3a","Metastatic")){
    sub_obj = subset(object, subset = stage == x)
    mat = sub_obj[["RNA"]]@data
    sub_result_list = list()
    for(x in 1:length(gene_list)){
      selected_genes = gene_list[[x]]
      marker_tag = names(gene_list[x])
      is_in_mat = selected_genes %in% rownames(mat)
      if(!all(is_in_mat)){
        not_available_genes = selected_genes[!is_in_mat]
        print(paste0(paste(not_available_genes, collapse = ","), 
                     " is not available for analysis"))
      }
      selected_genes = selected_genes[is_in_mat]
      selected_mat = mat[selected_genes, ]
      composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
      composite_marker_exp = scale(composite_marker_exp)
      marker_tag = paste0(marker_tag, z_score_tag)
      sub_result_list[[marker_tag]] = composite_marker_exp
    }
    sub_result_list = as.data.frame(sub_result_list)
    result_list = rbind(result_list, sub_result_list)
  }
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}
#add normalize count to object
AddNormalizedCount = function(object, gene_list){
  mat = object[["RNA"]]@data
  result_list = list()
  for(x in 1:length(gene_list)){
    selected_genes = gene_list[[x]]
    selected_genes = trimws(selected_genes, which = c("both"))
    marker_tag = names(gene_list[x])
    is_in_mat = selected_genes %in% rownames(mat)
    if(!all(is_in_mat)){
      not_available_genes = selected_genes[!is_in_mat]
      print(paste0(paste(not_available_genes, collapse = ","), 
                   " is not available for analysis"))
    }
    selected_genes = selected_genes[is_in_mat]
    selected_mat = mat[selected_genes, ]
    composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
    marker_tag = paste0(marker_tag,normalized_count_tag)
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}