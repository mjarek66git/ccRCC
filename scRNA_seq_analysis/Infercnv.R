prepare_infercnv_data = function(obj, annotation = "integrated_snn_res.0.1", 
                                 infercnv_dir = "./infercnv_inputs/"){
  
  # check if infercnv dir exists 
  if (!dir.exists(infercnv_dir)){
   dir.create(infercnv_dir, recursive = T) 
  }
  
  
  DefaultAssay(obj) = "RNA"
  #get count matrix
  counts_matrix = GetAssayData(obj, slot="counts")
  
  #get sample annotation
  sample_anno = obj@meta.data[, annotation]
  sample_anno_df = data.frame(cbind(colnames(counts_matrix), sample_anno))
  
  #get gene position anno
  human_position_anno = read.csv("human_gene_pos.csv", header = T)
  human_position_anno = human_position_anno[!duplicated(human_position_anno[c("gene_name")]),]
  all_gene = data.frame(cbind(rownames(counts_matrix), 1))
  all_gene = left_join(all_gene, human_position_anno, by=c("X1"="gene_name"))
  all_gene = na.omit(all_gene)
  all_gene = all_gene[,-2]
  
  #refine count matrix
  counts_matrix = counts_matrix[rownames(counts_matrix) %in% all_gene$X1, ]
  print(nrow(counts_matrix))
  
  #export gene and sample annotation
  sample_anno_df_path = paste0(infercnv_dir, "sample_annotation_df.tsv")
  gene_anno_df_path = paste0(infercnv_dir, "gene_annotation_df.tsv")
  write.table(sample_anno_df,sample_anno_df_path, sep = "\t", row.names = F, col.names = F)
  write.table(all_gene,gene_anno_df_path, sep = "\t", row.names = F, col.names = F)
  
  return(list(counts_matrix, sample_anno_df_path, gene_anno_df_path))
}