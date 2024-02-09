# # =======================================================================================
# =======================================================================================
# =======================================================================================
#541NO RL B1 = S1A  (Normal Tissue)

#541TU RL B2 = S1B (Tumor Tissue RL)
#543NO RL B7 = S1C (Normal Tissue)
#543NO TU B8 = S1D (Tumor Tissue RL)
#---------------------------------------------------
#257NO RL B9 = S2A (Normal Tissue)
#257TU RL B10 = S2B (Tumor Tissue RL)
#179TU DF B4 = S2C (Tumor Tissue DF)
#508TU DF B6 = S2D (Tumor Tissue DF)

# This analysis is provided in R version 4.3.2 (2023-10-31)

library(dplyr)
library(seewave)
library(circlize)
library(Seurat)
library(STutility)
library(MAST)
library(DOSE)
library(hash)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(kableExtra)
library(knitr)
library(stringr)
library(ComplexHeatmap)
library(pheatmap)
library(cowplot)
library("MALDIquant")
library(imagefx)
library(magick)
library(ggplot2)
library(ggcorrplot)
library(ggplot2)
library(ggridges)
library(corrplot)
library(hdf5r)
library(glmGamPoi)
library(stringr)
library(cluster)

# =========================================================================
# =========================================================================
# =========================================================================
# Part-1: Loading required functions
# =========================================================================
# =========================================================================
# =========================================================================
part_1 <- TRUE
if (part_1) {
  split_function <- function(inpute_item) {
    inpute_item_splitted <- str_split(inpute_item, "_")
    #print(inpute_item_splitted)
    #print(length(inpute_item_splitted[[1]]))
    clusttt <-
      inpute_item_splitted[[1]][length(inpute_item_splitted[[1]])]
    #print(clusttt)
    #print(substr(clusttt,1,nchar(clusttt) -1 ))
    #print("--------")
    return(substr(clusttt, 1, nchar(clusttt) - 1))
  }
  split_function("custom_features_scaled_top_20_01")
  
  # =======================================================================================
  # =======================================================================================
  # =======================================================================================
  
  AddCustomModuleScore <- function(object = object,
                                   features = features,
                                   name = name,
                                   ctrl = 100,
                                   nbin = 25,
                                   assay = assay,
                                   slot = "data",
                                   k = FALSE,
                                   seed = 1) {
    print(paste("calculating custom module score for", name))
    #print(assay)
    #print(slot)
    #pool = rownames(object@assays$assay@slot)
    # Find how many gene lists were provided. In this case just one.
    #cluster.length <- length(x = features)
    
    cluster.length <- 1
    
    # Pull the expression data from the provided Seurat object
    assay.data <-
      GetAssayData(object = object,
                   assay = assay,
                   slot = slot)
    
    #print(assay.data[1:5,1:5])
    # For all genes, get the average expression across all cells (named vector)
    data.avg <- Matrix::rowMeans(x = assay.data[,])
    # Order genes from lowest average expression to highest average expression
    data.avg <- data.avg[order(data.avg)]
    #print(data.avg)
    # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
    data.cut <-
      ggplot2::cut_number(
        x = data.avg + rnorm(n = length(data.avg)) / 1e+30,
        n = nbin,
        labels = FALSE,
        right = FALSE
      )
    
    # Set the names of the cuts as the gene names
    names(x = data.cut) <- names(x = data.avg)
    
    # Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
    ctrl.use <- c()
    
    # For each of the input gene lists:
    
    
    
    
    
    features.use_found <- c()
    features.use_not_found <- c()
    #cluster.length
    for (i in 1:cluster.length) {
      # Get the gene names from the input gene set as a character vector
      
      features.use <- unlist(features, recursive = FALSE)
      
      # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
      for (j in 1:length(x = features.use)) {
        # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
        #print("features.use[j]")
        #print(features.use[j])
        if (features.use[j] %in% names(data.cut)) {
          # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
          features.use_found <-
            append(features.use_found, features.use[j])
          sample_ctrl_list <-
            names(x = sample(
              x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
              size = ctrl,
              replace = FALSE
            ))
          for (sample_ctrl in sample_ctrl_list) {
            ctrl.use <- append(ctrl.use, sample_ctrl)
          }
          
          
        }
        else{
          features.use_not_found <-
            append(features.use_not_found, features.use[j])
          
        }
      }
    }
    print(paste(
      "features",
      features.use_not_found,
      "not found in the specified assay"
    ))
    #print("here1")
    
    # Have a quick look at what's in ctrl.use:
    #print(class(ctrl.use))
    ## [1] "list"
    #print(length(ctrl.use))
    ## [1] 1
    #print(class(ctrl.use))
    ## [1] "character"
    # There should be length(features.use)*ctrl genes (i.e. 20*100):
    #print(length(ctrl.use))
    #features.use <- features.use_found
    # =======================================================================
    # =======================================================================
    # =======================================================================
    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
    ctrl.use <- unlist(ctrl.use, recursive = FALSE)
    #print("here2")
    #print(length(ctrl.use))
    ## Get control gene scores
    
    # Create an empty matrix with dimensions;
    # number of rows equal to the number of gene sets (just one here)
    # number of columns equal to number of cells in input Seurat object
    ctrl.scores <- matrix(
      data = numeric(length = 1L),
      nrow = length(x = ctrl.use),
      ncol = ncol(x = object)
    )
    #print("here3")
    #print(dim(ctrl.scores))
    #print(ctrl.use)
    # Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
    #for (i in 1:length(ctrl.use)) {
    # Get control gene names as a vector
    features.use <- ctrl.use
    # For each cell, calculate the mean expression of *all* of the control genes
    #print(dim(assay.data[features.use,]))
    ctrl.scores <- Matrix::colMeans(x = assay.data[features.use, ])
    #}
    
    #print("here4")
    #print(ctrl.scores)
    ## Get scores for input gene sets
    
    # Similar to the above, create an empty matrix
    features.scores <- matrix(
      data = numeric(length = 1L),
      nrow = cluster.length,
      ncol = ncol(x = object)
    )
    #print("here5")
    #print(cluster.length)
    #print(features.scores)
    # Loop through input gene sets and calculate the mean expression of these genes for each cell
    for (i in 1:cluster.length) {
      features.use <- unlist(features.use_found, recursive = FALSE)
      #print(features.use)
      data.use <- assay.data[features.use, ]
      #print(data.use)
      features.scores[i,] <- Matrix::colMeans(x = data.use)
      #print(features.scores[i, ] )
    }
    #print("here6")
    # =======================================================================
    # =======================================================================
    # =======================================================================
    #print(dim(features.scores))
    #print(dim(ctrl.scores))
    features.scores.use <- features.scores - ctrl.scores
    #print("here7")
    # Name the result the "name" variable + whatever the position the geneset was in the input list, e.g. "Cluster1"
    rownames(x = features.scores.use) <-
      paste0(name, 1:cluster.length)
    #print("here8")
    # Change the matrix from wide to long
    features.scores.use <-
      as.data.frame(x = t(x = features.scores.use))
    #print("here9")
    # Give the rows of the matrix, the names of the cells
    rownames(x = features.scores.use) <- colnames(x = object)
    #print("here10")
    #print(features.scores.use)
    # Add the result as a metadata column to the input Seurat object
    object[[colnames(x = features.scores.use)]] <- features.scores.use
    #print("here11")
    #print(colnames(x = features.scores.use))
    return(features.scores.use)
    # Voila!
    #FeaturePlot(object,
    #            features = "custom_tca_features1") +
    #  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
    #print("here12")
  }
  
  
  # ========================================================================
  # ========================================================================
  # ========================================================================
  
  split_function_immune_state <- function(inpute_item) {
    inpute_item_splitted <- str_split(inpute_item, "_")
    #print(inpute_item_splitted)
    #print(length(inpute_item_splitted[[1]]))
    clusttt <-
      inpute_item_splitted[[1]][length(inpute_item_splitted[[1]]) - 1]
    #print(clusttt)
    #print(substr(clusttt,1,nchar(clusttt) -1 ))
    #print("--------")
    return(clusttt)
  }
  
  
  collocalization_post_process_four_sections <-
    function(object,
             neighbors_in_object,
             item_name1,
             item_name2,
             meta_item) {
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      my_list <- rep(0.0, nrow(object[[]]))
      
      #for(row in 1:50) {
      for (row in 1:nrow(neighbors_in_object)) {
        #print(row)
        cal_value1 <- 0.0
        cal_value2 <- 0.0
        cell_id <- row.names(neighbors_in_object)[row]
        
        if (cell_id %in% rr_nn) {
          if ((object[[]][cell_id, item_name1] * object[[]][cell_id, item_name2]) < 0.0) {
            print("here1")
            flag1 = TRUE
            flag2 = FALSE
            if (object[[]][cell_id, item_name1] > object[[]][cell_id, item_name2]) {
              cal_value1 <-
                cal_value2 + abs(object[[]][cell_id, item_name1]) + abs(object[[]][cell_id, item_name2])
            }
            else{
              cal_value1 <-
                cal_value2 - abs(object[[]][cell_id, item_name1]) - abs(object[[]][cell_id, item_name2])
            }
          }
          else{
            flag2 = TRUE
            flag1 = FALSE
            #print("else1")
            if (object[[]][cell_id, item_name1] > 0.0) {
              cal_value2 <-
                cal_value2  + abs(object[[]][cell_id, item_name1]) + abs(object[[]][cell_id, item_name2])
            }
            else{
              cal_value2 <-
                cal_value2  - abs(object[[]][cell_id, item_name1]) - abs(object[[]][cell_id, item_name2])
            }
          }
          
          
          
          
          
          #if(object[[]][cell_id,item_name] > threashold)
          #{
          iter <- 0
          for (col in 1:ncol(neighbors_in_object)) {
            cell_id2 <-
              row.names(neighbors_in_object)[neighbors_in_object[row, col]]
            #print("col")
            #print(row.names(neighbor_df_2D)[col])
            #print(row.names(neighbor_df_2D)[col])
            if (cell_id2 %in% rr_nn) {
              if ((object[[]][cell_id2, item_name1] * object[[]][cell_id2, item_name2]) < 0.0) {
                print("here2")
                if (object[[]][cell_id2, item_name1] > object[[]][cell_id2, item_name2]) {
                  cal_value2 <-
                    cal_value2 + abs(object[[]][cell_id2, item_name1]) + abs(object[[]][cell_id2, item_name2])
                }
                else{
                  cal_value2 <-
                    cal_value2 - abs(object[[]][cell_id2, item_name1]) - abs(object[[]][cell_id2, item_name2])
                }
              }
              else{
                #print("else2")
                if (object[[]][cell_id2, item_name1] > 0.0) {
                  cal_value2 <-
                    cal_value2 + abs(object[[]][cell_id2, item_name1]) + abs(object[[]][cell_id2, item_name2])
                }
                else{
                  cal_value2 <-
                    cal_value2 - abs(object[[]][cell_id2, item_name1]) - abs(object[[]][cell_id2, item_name2])
                }
                
              }
              
              
            }
            
          }
          
          p1 <- match(cell_id, rr_nn)
          #print("p1")
          #print(p1)
          
          my_list[p1] <- cal_value2 / 14.0
          
          
          
        }
      }
      #print(my_list)
      #print(meta_item)
      CellsMeta[meta_item] <- my_list
      
      
      object <- AddMetaData(object, CellsMeta)
      return(object)
    }
  
  
  neighborhood_histogram_analysis_for_one_cluster_reverse <-
    function(object,
             neighbors_in_object,
             input_cluster_id) {
      #object <- tumor_1B_subset
      #neighbors_in_object <- neighbor_df_1B
      #input_cluster_id <- 0
      #head(object[[]])
      #object[[]]["TATGGCAGACTTTCGA-1_1",]
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      col_class_vector <- c()
      col_value_vector <- c()
      col_class_vector_decoy <- c()
      col_value_vector_decoy <- c()
      decoy_vec <- c()
      decoy_vec2 <- c()
      error_iter <- 0
      cluster_list <- unique(object[[]][, "seurat_clusters"])
      first_row <- TRUE
      
      
      #new_row = c(id = 33, pages=50, name = "java")
      #df_orig = rbind(df,new_row)
      #new_row = c(id = 33, pages=50, name = "java")
      #df_decoy = rbind(df,new_row)
      
      
      for (iter in 1:length(summary(object[[]][, "seurat_clusters"]))) {
        classt <- names(summary(object[[]][, "seurat_clusters"]))[iter]
        rept <- summary(object[[]][, "seurat_clusters"])[[iter]]
        if (classt != input_cluster_id) {
          #print(iter)
          #print(class)
          #print(rept)
          #print("----")
          decoy_vec2 <- append(decoy_vec2, rep(classt, rept))
        }
        
      }
      #table(decoy_vec)
      x_rand2 <- sample(decoy_vec2)        # Sample vector
      decoy_vec2_iter <- 1
      
      for (iter in 1:length(object[[]][, "seurat_clusters"])) {
        classt <- object[[]][, "seurat_clusters"][iter]
        
        if (classt != input_cluster_id) {
          #print(iter)
          #print(class)
          #print(rept)
          #print("----")
          decoy_vec <-
            append(decoy_vec, as.integer(x_rand2[decoy_vec2_iter]))
          decoy_vec2_iter <- decoy_vec2_iter + 1
        } else{
          decoy_vec <- append(decoy_vec, input_cluster_id)
        }
        
      }
      
      print(table(decoy_vec))
      dim(object[[]])
      empty_orig_Matrix <-
        matrix(0, ncol = length(unique(object[[]][, "seurat_clusters"])), nrow = 1)
      empty_decoy_Matrix <-
        matrix(0, ncol = length(unique(object[[]][, "seurat_clusters"])), nrow = 1)
      col_names_orig <- c()
      col_names_decoy <- c()
      for (cl in cluster_list) {
        print(cl)
        Newcolname <- paste0("orig_", cl)
        col_names_orig <- append(col_names_orig, Newcolname)
        #df_orig_row[[Newcolname]] <- 0
        Newcolname <- paste0("decoy_", cl)
        col_names_decoy <- append(col_names_decoy, Newcolname)
        #df_decoy_row[[Newcolname]] <- 0
      }
      colnames(empty_orig_Matrix) <- col_names_orig
      colnames(empty_decoy_Matrix) <- col_names_decoy
      
      print(table(object[[]][, "seurat_clusters"]))
      print(table(object[[]][, "seurat_clusters"]))
      for (row in 1:nrow(neighbors_in_object)) {
        if (row %% 500 == 0) {
          print(row)
        }
        
        cell_id <- row.names(neighbors_in_object)[row]
        
        column <- object[[]][cell_id, "seurat_clusters"]
        column_name <- paste0(column)
        
        
        if (cell_id %in% rr_nn) {
          row_decoy <- which(cell_id == rr_nn)[[1]]
          
          column_decoy <- as.numeric(decoy_vec[row_decoy])
          
          
          column_name_decoy <- paste0(column_decoy)
          
          
          
          df_orig_row <- data.frame(empty_orig_Matrix)
          df_decoy_row <- data.frame(empty_decoy_Matrix)
          
          #print(df_decoy_row)
          #print(df_orig_row)
          iter <- 1.0
          iter_decoy <- 1.0
          print("column")
          print(column)
          print("column_decoy")
          print(column_decoy)
          print("++++++++")
          if (column != input_cluster_id) {
            print("column")
            print(column)
            print("------")
            for (col in 1:ncol(neighbors_in_object)) {
              cell_id2 <-
                row.names(neighbors_in_object)[neighbors_in_object[row, col]]
              #print("col")
              #print(row.names(neighbor_df_2D)[col])
              #print(row.names(neighbor_df_2D)[col])
              if (cell_id2 %in% rr_nn) {
                iter <- iter + 1.0
                
                nei_cl_orig = object[[]][cell_id2, "seurat_clusters"]
                
                if (nei_cl_orig == input_cluster_id) {
                  Newcolname <- paste0("orig_", column)
                  df_orig_row[[Newcolname]] <-
                    df_orig_row[Newcolname] + 1
                  
                }
              }
              #}else{
              # error_iter <- error_iter + 1
              #}
              
            }
          }
          
          if (column_decoy != input_cluster_id) {
            print("column_decoy")
            print(column_decoy)
            print("------")
            for (col in 1:ncol(neighbors_in_object)) {
              cell_id2 <-
                row.names(neighbors_in_object)[neighbors_in_object[row, col]]
              #print("col")
              #print(row.names(neighbor_df_2D)[col])
              #print(row.names(neighbor_df_2D)[col])
              if (cell_id2 %in% rr_nn) {
                row_decoy2 <- which(cell_id2 == rr_nn)[[1]]
                
                iter_decoy <- iter_decoy + 1.0
                
                nei_cl_decoy = decoy_vec[row_decoy2]
                
                if (nei_cl_decoy == input_cluster_id) {
                  Newcolname <- paste0("decoy_", column_decoy)
                  df_decoy_row[[Newcolname]] <-
                    df_decoy_row[Newcolname] + 1
                  
                }
              }
              #}else{
              # error_iter <- error_iter + 1
              #}
              
            }
          }
          
          if (first_row) {
            df_orig <- df_orig_row
            df_decoy <- df_decoy_row
            first_row <- FALSE
          } else{
            df_orig = rbind(df_orig, df_orig_row)
            df_decoy = rbind(df_decoy, df_decoy_row)
            
          }
          
          #new_row = c(id = 33, pages=50, name = "java")
          #df_orig = rbind(df,new_row)
          
          #col_value_vector <- append(col_value_vector, cal_value/iter)
          #col_value_vector_decoy <- append(col_value_vector_decoy, cal_value_decoy/iter_decoy)
          
          
        }
      }
      
      
      
      
      #print(dim(df_orig))
      #print(dim(df_decoy))
      updated <- data.frame()
      updated <-
        data.frame(
          "orig_0" = as.numeric(unlist(df_orig$orig_0)),
          "orig_1" = as.numeric(unlist(df_orig$orig_1)),
          "orig_2" = as.numeric(unlist(df_orig$orig_2)),
          "orig_3" = as.numeric(unlist(df_orig$orig_3)),
          "orig_4" = as.numeric(unlist(df_orig$orig_4)),
          "orig_5" = as.numeric(unlist(df_orig$orig_5)),
          "orig_6" = as.numeric(unlist(df_orig$orig_6)),
          "orig_7" = as.numeric(unlist(df_orig$orig_7)),
          "orig_8" = as.numeric(unlist(df_orig$orig_8)),
          "decoy_0" = as.numeric(unlist(df_decoy$decoy_0)),
          "decoy_1" = as.numeric(unlist(df_decoy$decoy_1)),
          "decoy_2" = as.numeric(unlist(df_decoy$decoy_2)),
          "decoy_3" = as.numeric(unlist(df_decoy$decoy_3)),
          "decoy_4" = as.numeric(unlist(df_decoy$decoy_4)),
          "decoy_5" = as.numeric(unlist(df_decoy$decoy_5)),
          "decoy_6" = as.numeric(unlist(df_decoy$decoy_6)),
          "decoy_7" = as.numeric(unlist(df_decoy$decoy_7)),
          "decoy_8" = as.numeric(unlist(df_decoy$decoy_8))
        )
      
      #updated <- cbind(df_orig, df_decoy)
      if (FALSE) {
        my_list <- list(class = col_class_vector,
                        class_decoy = decoy_vec)
        updated <- cbind(updated, my_list)
      }
      
      
      #return(data.frame(updated))
      return(updated)
    }
  
  
  
  neighborhood_histogram_analysis_for_one_cluster <-
    function(object,
             neighbors_in_object,
             input_cluster_id) {
      #object <- tumor_1B_subset
      #neighbors_in_object <- neighbor_df_1B
      #input_cluster_id <- 0
      #head(object[[]])
      #object[[]]["TATGGCAGACTTTCGA-1_1",]
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      col_class_vector <- c()
      col_value_vector <- c()
      col_class_vector_decoy <- c()
      col_value_vector_decoy <- c()
      decoy_vec <- c()
      decoy_vec2 <- c()
      error_iter <- 0
      cluster_list <- unique(object[[]][, "seurat_clusters"])
      first_row <- TRUE
      
      
      #new_row = c(id = 33, pages=50, name = "java")
      #df_orig = rbind(df,new_row)
      #new_row = c(id = 33, pages=50, name = "java")
      #df_decoy = rbind(df,new_row)
      
      
      for (iter in 1:length(summary(object[[]][, "seurat_clusters"]))) {
        classt <- names(summary(object[[]][, "seurat_clusters"]))[iter]
        rept <- summary(object[[]][, "seurat_clusters"])[[iter]]
        if (classt != input_cluster_id) {
          #print(iter)
          #print(class)
          #print(rept)
          #print("----")
          decoy_vec2 <- append(decoy_vec2, rep(classt, rept))
        }
        
      }
      #table(decoy_vec)
      x_rand2 <- sample(decoy_vec2)        # Sample vector
      decoy_vec2_iter <- 1
      
      for (iter in 1:length(object[[]][, "seurat_clusters"])) {
        classt <- object[[]][, "seurat_clusters"][iter]
        
        if (classt != input_cluster_id) {
          #print(iter)
          #print(class)
          #print(rept)
          #print("----")
          decoy_vec <- append(decoy_vec, x_rand2[decoy_vec2_iter])
          decoy_vec2_iter <- decoy_vec2_iter + 1
        } else{
          decoy_vec <- append(decoy_vec, input_cluster_id)
        }
        
      }
      
      print(table(decoy_vec))
      dim(object[[]])
      empty_orig_Matrix <-
        matrix(0, ncol = length(unique(object[[]][, "seurat_clusters"])), nrow = 1)
      empty_decoy_Matrix <-
        matrix(0, ncol = length(unique(object[[]][, "seurat_clusters"])), nrow = 1)
      col_names_orig <- c()
      col_names_decoy <- c()
      for (cl in cluster_list) {
        print(cl)
        Newcolname <- paste0("orig_", cl)
        col_names_orig <- append(col_names_orig, Newcolname)
        #df_orig_row[[Newcolname]] <- 0
        Newcolname <- paste0("decoy_", cl)
        col_names_decoy <- append(col_names_decoy, Newcolname)
        #df_decoy_row[[Newcolname]] <- 0
      }
      colnames(empty_orig_Matrix) <- col_names_orig
      colnames(empty_decoy_Matrix) <- col_names_decoy
      
      
      for (row in 1:nrow(neighbors_in_object)) {
        if (row %% 500 == 0) {
          print(row)
        }
        
        cell_id <- row.names(neighbors_in_object)[row]
        
        column <- object[[]][cell_id, "seurat_clusters"]
        column_name <- paste0(column)
        
        
        if (cell_id %in% rr_nn) {
          if (column == input_cluster_id) {
            row_decoy <- which(cell_id == rr_nn)[[1]]
            column_decoy <- decoy_vec[row_decoy]
            
            column_name_decoy <- paste0(column_decoy)
            
            col_class_vector <- append(col_class_vector, column_name)
            
            col_class_vector_decoy <-
              append(col_class_vector_decoy, column_name_decoy)
            
            df_orig_row <- data.frame(empty_orig_Matrix)
            df_decoy_row <- data.frame(empty_decoy_Matrix)
            
            #print(df_decoy_row)
            #print(df_orig_row)
            iter <- 1.0
            iter_decoy <- 1.0
            
            for (col in 1:ncol(neighbors_in_object)) {
              cell_id2 <-
                row.names(neighbors_in_object)[neighbors_in_object[row, col]]
              #print("col")
              #print(row.names(neighbor_df_2D)[col])
              #print(row.names(neighbor_df_2D)[col])
              if (cell_id2 %in% rr_nn) {
                iter <- iter + 1.0
                
                nei_cl_orig = object[[]][cell_id2, "seurat_clusters"]
                
                
                row_decoy2 <- which(cell_id2 == rr_nn)[[1]]
                
                iter_decoy <- iter_decoy + 1.0
                
                nei_cl_decoy = decoy_vec[row_decoy2]
                
                Newcolname <- paste0("orig_", nei_cl_orig)
                df_orig_row[[Newcolname]] <- df_orig_row[Newcolname] + 1
                Newcolname <- paste0("decoy_", nei_cl_decoy)
                df_decoy_row[[Newcolname]] <- df_decoy_row[Newcolname] + 1
                #print(df_orig_row)
                #print(df_decoy_row)
                #print("-----")
              }
              #}else{
              # error_iter <- error_iter + 1
              #}
              
            }
            
            if (first_row) {
              df_orig <- df_orig_row
              df_decoy <- df_decoy_row
              first_row <- FALSE
            } else{
              df_orig = rbind(df_orig, df_orig_row)
              df_decoy = rbind(df_decoy, df_decoy_row)
              
            }
            
            #new_row = c(id = 33, pages=50, name = "java")
            #df_orig = rbind(df,new_row)
            
            #col_value_vector <- append(col_value_vector, cal_value/iter)
            #col_value_vector_decoy <- append(col_value_vector_decoy, cal_value_decoy/iter_decoy)
            
          }
        }
      }
      
      
      
      
      #print(dim(df_orig))
      #print(dim(df_decoy))
      updated <- data.frame()
      updated <-
        data.frame(
          "orig_0" = as.numeric(unlist(df_orig$orig_0)),
          "orig_1" = as.numeric(unlist(df_orig$orig_1)),
          "orig_2" = as.numeric(unlist(df_orig$orig_2)),
          "orig_3" = as.numeric(unlist(df_orig$orig_3)),
          "orig_4" = as.numeric(unlist(df_orig$orig_4)),
          "orig_5" = as.numeric(unlist(df_orig$orig_5)),
          "orig_6" = as.numeric(unlist(df_orig$orig_6)),
          "orig_7" = as.numeric(unlist(df_orig$orig_7)),
          "orig_8" = as.numeric(unlist(df_orig$orig_8)),
          "decoy_0" = as.numeric(unlist(df_decoy$decoy_0)),
          "decoy_1" = as.numeric(unlist(df_decoy$decoy_1)),
          "decoy_2" = as.numeric(unlist(df_decoy$decoy_2)),
          "decoy_3" = as.numeric(unlist(df_decoy$decoy_3)),
          "decoy_4" = as.numeric(unlist(df_decoy$decoy_4)),
          "decoy_5" = as.numeric(unlist(df_decoy$decoy_5)),
          "decoy_6" = as.numeric(unlist(df_decoy$decoy_6)),
          "decoy_7" = as.numeric(unlist(df_decoy$decoy_7)),
          "decoy_8" = as.numeric(unlist(df_decoy$decoy_8))
        )
      
      #updated <- cbind(df_orig, df_decoy)
      if (FALSE) {
        my_list <- list(class = col_class_vector,
                        class_decoy = decoy_vec)
        updated <- cbind(updated, my_list)
      }
      
      
      #return(data.frame(updated))
      return(updated)
    }
  
  
  neighborhood_histogram_analysis <-
    function(object, neighbors_in_object) {
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      col_class_vector <- c()
      col_value_vector <- c()
      col_class_vector_decoy <- c()
      col_value_vector_decoy <- c()
      decoy_vec <- c()
      error_iter <- 0
      
      for (iter in 1:length(summary(object[[]][, "seurat_clusters"]))) {
        classt <- names(summary(object[[]][, "seurat_clusters"]))[iter]
        rept <- summary(object[[]][, "seurat_clusters"])[[iter]]
        print(iter)
        print(class)
        print(rept)
        print("----")
        decoy_vec <- append(decoy_vec, rep(classt, rept))
        
      }
      
      #table(decoy_vec)
      x_rand <- sample(decoy_vec)        # Sample vector
      
      for (row in 1:nrow(neighbors_in_object)) {
        cell_id <- row.names(neighbors_in_object)[row]
        column <- object[[]][cell_id, "seurat_clusters"]
        column_name <- paste0("sp_", column)
        
        
        
        if (cell_id %in% rr_nn) {
          row_decoy <- which(cell_id == rr_nn)[[1]]
          column_decoy <- x_rand[row_decoy]
          
          column_name_decoy <- paste0("sp_", column_decoy)
          
          col_class_vector <- append(col_class_vector, column_name)
          
          col_class_vector_decoy <-
            append(col_class_vector_decoy, column_name_decoy)
          
          cal_value <- 1.0
          iter <- 1.0
          
          cal_value_decoy <- 1.0
          iter_decoy <- 1.0
          
          for (col in 1:ncol(neighbors_in_object)) {
            cell_id2 <-
              row.names(neighbors_in_object)[neighbors_in_object[row, col]]
            #print("col")
            #print(row.names(neighbor_df_2D)[col])
            #print(row.names(neighbor_df_2D)[col])
            if (cell_id2 %in% rr_nn) {
              iter <- iter + 1.0
              
              if (column == object[[]][cell_id2, "seurat_clusters"]) {
                cal_value <- cal_value + 1.0
                
              }
              #}
              #if ( neighbors_in_object[row_decoy,col_decoy] < length(x_rand) + 1){
              row_decoy2 <- which(cell_id2 == rr_nn)[[1]]
              #print(neighbors_in_object[row_decoy,col])
              #print(column_decoy)
              #print(x_rand[neighbors_in_object[row_decoy,col]])
              #print("-----")
              iter_decoy <- iter_decoy + 1.0
              
              if (column_decoy == x_rand[row_decoy2]) {
                cal_value_decoy <- cal_value_decoy + 1.0
                
              }
            }
            #}else{
            # error_iter <- error_iter + 1
            #}
            
          }
          col_value_vector <- append(col_value_vector, cal_value / iter)
          col_value_vector_decoy <-
            append(col_value_vector_decoy, cal_value_decoy / iter_decoy)
          
        }
      }
      
      my_list <- list(
        class = col_class_vector,
        value = col_value_vector,
        class_decoy = col_class_vector_decoy,
        value_decoy = col_value_vector_decoy
      )
      print("error_iter:")
      print(error_iter)
      return(data.frame(my_list))
    }
  
  
  collocalization_post_process_summed <-
    function(object,
             neighbors_in_object,
             item_name1,
             item_name2,
             meta_item) {
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      my_list <- rep(0.0, nrow(object[[]]))
      my_list2 <- rep(0.0, nrow(object[[]]))
      
      #for(row in 1:50) {
      for (row in 1:nrow(neighbors_in_object)) {
        #print(row)
        cal_value <- 0.0
        cal_value2 <- 0.0
        cell_id <- row.names(neighbors_in_object)[row]
        
        if (cell_id %in% rr_nn) {
          #cal_value <- cal_value + object[[]][cell_id,item_name1] + object[[]][cell_id,item_name2]
          
          
          if ((object[[]][cell_id, item_name1] * object[[]][cell_id, item_name2]) < 0.0) {
            dis = TRUE
            con = FALSE
            if (object[[]][cell_id, item_name1] > object[[]][cell_id, item_name2]) {
              cal_value2 <-
                cal_value2 + abs(object[[]][cell_id, item_name1]) + abs(object[[]][cell_id, item_name2])
            }
            else{
              cal_value2 <-
                cal_value2 - abs(object[[]][cell_id, item_name1]) - abs(object[[]][cell_id, item_name2])
            }
          }
          else{
            con = TRUE
            dis = FALSE
            if (object[[]][cell_id, item_name1] > 0.0) {
              cal_value <-
                cal_value + abs(object[[]][cell_id, item_name1]) + abs(object[[]][cell_id, item_name2])
            }
            else{
              cal_value <-
                cal_value - abs(object[[]][cell_id, item_name1]) - abs(object[[]][cell_id, item_name2])
            }
            
          }
          
          
          
          
          #if(object[[]][cell_id,item_name] > threashold)
          #{
          iter <- 0
          for (col in 1:ncol(neighbors_in_object)) {
            cell_id2 <-
              row.names(neighbors_in_object)[neighbors_in_object[row, col]]
            #print("col")
            #print(row.names(neighbor_df_2D)[col])
            #print(row.names(neighbor_df_2D)[col])
            if (cell_id2 %in% rr_nn) {
              #cal_value <- cal_value + object[[]][cell_id2,item_name1] + object[[]][cell_id2,item_name2]
              
              if ((object[[]][cell_id2, item_name1] * object[[]][cell_id2, item_name2]) < 0.0) {
                if (object[[]][cell_id2, item_name1] > object[[]][cell_id2, item_name2]) {
                  cal_value2 <-
                    cal_value2 + abs(object[[]][cell_id2, item_name1]) + abs(object[[]][cell_id2, item_name2])
                }
                else{
                  cal_value2 <-
                    cal_value2 - abs(object[[]][cell_id2, item_name1]) - abs(object[[]][cell_id2, item_name2])
                }
              }
              else{
                if (object[[]][cell_id2, item_name1] > 0.0) {
                  cal_value <-
                    cal_value + abs(object[[]][cell_id2, item_name1]) + abs(object[[]][cell_id2, item_name2])
                }
                else{
                  cal_value <-
                    cal_value - abs(object[[]][cell_id2, item_name1]) - abs(object[[]][cell_id2, item_name2])
                }
                
              }
              
              
            }
            
          }
          
          p1 <- match(cell_id, rr_nn)
          #print("p1")
          #print(p1)
          if (con) {
            my_list[p1] <- cal_value / 14.0
            my_list2[p1] <- 0.0
          }
          else{
            my_list[p1] <- 0.0
            my_list2[p1] <- cal_value2 / 14.0
            
          }
          
          
          
        }
      }
      #print(my_list)
      #print(meta_item)
      CellsMeta[meta_item] <- my_list
      
      discordant <- paste0("dis_", meta_item)
      CellsMeta[discordant] <- my_list2
      
      object <- AddMetaData(object, CellsMeta)
      return(object)
    }
  
  
  
  spatial_cluster_neighborhood_analysis <-
    function(object,
             neighbors_in_object,
             metadata_col,
             cluster_number,
             nneighbor) {
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      my_list <- rep(0.0, nrow(object[[]]))
      my_list2 <- rep(0.0, nrow(object[[]]))
      #object <-tumor_1D_subset
      cluster_hash <- hash()
      unique_spots <- c()
      #metadata_col <- "seurat_clusters"
      for (name_c in names(table(object[[]][, metadata_col]))) {
        cluster_hash[name_c] <- 0
      }
      
      for (row in 1:nrow(neighbors_in_object)) {
        cell_id <- row.names(neighbors_in_object)[row]
        
        if (cell_id %in% rr_nn &
            object[[]][cell_id, metadata_col] == cluster_number) {
          counter = 0
          for (col in 1:ncol(neighbors_in_object)) {
            cell_id2 <-
              row.names(neighbors_in_object)[neighbors_in_object[row, col]]
            
            if (cell_id2 %in% rr_nn &
                object[[]][cell_id2, metadata_col] == cluster_number) {
              counter <- counter + 1
            }
          }
          if (counter > nneighbor) {
            for (col in 1:ncol(neighbors_in_object)) {
              cell_id2 <-
                row.names(neighbors_in_object)[neighbors_in_object[row, col]]
              
              if (cell_id2 %in% rr_nn & !(cell_id2 %in% unique_spots)) {
                unique_spots <- append(unique_spots, cell_id2)
                cluster_class <-
                  toString(object[[]][cell_id2, metadata_col])
                
                cluster_hash[cluster_class] <-
                  cluster_hash[[cluster_class]] + 1
                
              }
            }
          }
          
        }
      }
      
      
      return(cluster_hash)
    }
  
  
  collocalization_post_process <-
    function(object,
             neighbors_in_object,
             item_name,
             threashold) {
      rr_nn <- row.names(object[[]])
      CellsMeta = object[[]]
      #print("here")
      my_list <- rep(0, nrow(object[[]]))
      #for(row in 1:50) {
      for (row in 1:nrow(neighbors_in_object)) {
        #print(row)
        
        cell_id <- row.names(neighbors_in_object)[row]
        if (cell_id %in% rr_nn) {
          if (object[[]][cell_id, item_name] > threashold)
          {
            iter <- 0
            for (col in 1:ncol(neighbors_in_object)) {
              cell_id2 <-
                row.names(neighbors_in_object)[neighbors_in_object[row, col]]
              #print("col")
              #print(row.names(neighbor_df_2D)[col])
              #print(row.names(neighbor_df_2D)[col])
              
              if (cell_id2 %in% rr_nn) {
                if (object[[]][cell_id2, item_name] > threashold)
                {
                  iter <- iter + 1
                }
                
              }
              
            }
            if (iter > 2) {
              #print(cell_id)
              #print(rr_nn)
              p1 <- match(cell_id, rr_nn)
              #print("p1")
              #print(p1)
              my_list[p1] <- object[[]][cell_id, item_name]
            }
            else{
              p1 <- match(cell_id, rr_nn)
              my_list[p1] <- -1.0 * abs(object[[]][cell_id, item_name])
            }
            
            
          }
          else{
            #print(cell_id)
            #print(rr_nn)
            p1 <- match(cell_id, rr_nn)
            #print("p1")
            #print(p1)
            my_list[p1] <- -1.0 * abs(object[[]][cell_id, item_name])
          }
          
          
        }
        
        
      }
      CellsMeta[paste0("colloc_", item_name)] <- my_list
      object <- AddMetaData(object, CellsMeta)
      return(object)
    }
  
  find_neighbor_of_spots <-
    function(object) {
      # create a function with the name my_function
      
      st.object <- object@tools$Staffli
      
      # Check if images are loaded
      if (length(x = st.object@rasterlists) == 0)
        stop("The images need to be loaded before running this function ... \n")
      
      # spatial information
      xys = setNames(st.object@meta.data[, c("pixel_x", "pixel_y", "sample")], c("x", "y", "sample"))
      unique(xys$sample)
      # Split x, y, s
      xys.list <- split(xys, xys$sample)[unique(xys$sample)]
      head(xys.list)
      # Obtain platforms
      platforms <- st.object@platforms
      platforms
      nNeighbours <- 6
      maxdist <- 1.5
      i = 1
      xys <- xys.list[[i]]
      print(i)
      print(xys)
      # vector matching spot_ID and order
      spotnames <- rownames(xys)
      print(spotnames)
      names(spotnames) <- c(1:nrow(xys)) %>% paste0()
      print(names(spotnames))
      # Get spot distances
      sdist <- st.object@pixels.per.um[i]
      sdist
      nNeighbours <-
        nNeighbours %||% ifelse(platforms[i] == "Visium", 6, 4)
      maxdist <-
        maxdist %||% ifelse(platforms[i] == "1k", 270 * sdist, 150 * sdist)
      nNeighbours
      maxdist
      if (!requireNamespace("dbscan"))
        install.packages("dbscan")
      knn_spatial <-
        dbscan::kNN(x = xys[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
      return(knn_spatial$id)
      
    }
  
  
  color_schema <- readRDS("~/Documents/color_schema.rds")
  color_schema
  options(future.globals.maxSize = 4000 * 1024 ^ 5)
  if_print = FALSE
  spatial_states_hash <- hash()
  
  
  # ---------------------------------------------------------------------------------------
  # immune cell gene markers
  # ---------------------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------------------
  super_genes <-
    read.csv("~/Documents/SuperGeneSets_10_2022.csv",
             header = T,
             sep = ",")
  colnames(super_genes)
  super_genes[["glycolysis"]]
  super_genes[["HIF1A"]]
  super_genes_hash <- hash()
  super_genes_hash_complex <- hash()
  super_genes_hash_complex[["complex.I"]] <-
    super_genes$complex.I[super_genes$complex.I != ""]
  super_genes_hash_complex[["complex.II"]] <-
    super_genes$complex.II[super_genes$complex.II != ""]
  super_genes_hash_complex[["complex.III"]] <-
    super_genes$complex.III[super_genes$complex.III != ""]
  super_genes_hash_complex[["complex.IV"]] <-
    super_genes$complex.IV[super_genes$complex.IV != ""]
  #super_genes_hash_complex[["Complex.V"]] <- super_genes$Complex.V_v1[super_genes$Complex.V_v1 != ""]
  super_genes_hash_complex[["Complex.V"]] <-
    super_genes$Complex.V_v2[super_genes$Complex.V_v2 != ""]
  super_genes_hash_complex
  oxPhos_gs <- super_genes$complex.I[super_genes$complex.I != ""]
  
  for (item in super_genes$complex.II[super_genes$complex.II != ""]) {
    oxPhos_gs <- append(oxPhos_gs, item)
  }
  for (item in super_genes$complex.III[super_genes$complex.III != ""]) {
    oxPhos_gs <- append(oxPhos_gs, item)
  }
  for (item in super_genes$complex.IV[super_genes$complex.IV != ""]) {
    oxPhos_gs <- append(oxPhos_gs, item)
  }
  for (item in super_genes$Complex.V_v2[super_genes$Complex.V_v2 != ""]) {
    oxPhos_gs <- append(oxPhos_gs, item)
  }
  names(super_genes)
  super_genes$complex.I
  super_genes_hash[["oxPhos"]] <- oxPhos_gs
  super_genes_hash[["TCA"]] <- super_genes$TCA[super_genes$TCA != ""]
  super_genes_hash[["glycolysis"]] <-
    super_genes$glycolysis[super_genes$glycolysis != ""]
  super_genes_hash[["MAS"]] <- super_genes$MAS[super_genes$MAS != ""]
  super_genes_hash[["HIF1A"]] <-
    super_genes$HIF1A[super_genes$HIF1A != ""]
  super_genes_hash[["Cu"]] <- super_genes$Cu[super_genes$Cu != ""]
  super_genes_hash[["EMT"]] <- super_genes$EMT[super_genes$EMT != ""]
  super_genes_hash[["HIF1A.2A"]] <-
    super_genes$HIF1A.2A[super_genes$HIF1A.2A != ""]
  super_genes_hash[["MTs"]] <- super_genes$MTs[super_genes$MTs != ""]
  super_genes_hash[["FAO"]] <- super_genes$FAO[super_genes$FAO != ""]
  super_genes_hash[["FAS"]] <- super_genes$FAS[super_genes$FAS != ""]
  #super_genes_hash[["Plamalogens"]] <- super_genes$Plamalogens[super_genes$Plamalogens != ""]
  super_genes_hash[["NRF2_targets"]] <-
    super_genes$NRF2[super_genes$NRF2 != ""]
  super_genes_hash[["complex.I"]] <-
    super_genes$complex.I[super_genes$complex.I != ""]
  super_genes_hash[["complex.II"]] <-
    super_genes$complex.II[super_genes$complex.II != ""]
  super_genes_hash[["complex.III"]] <-
    super_genes$complex.III[super_genes$complex.III != ""]
  super_genes_hash[["complex.IV"]] <-
    super_genes$complex.IV[super_genes$complex.IV != ""]
  #super_genes_hash_complex[["Complex.V"]] <- super_genes$Complex.V_v1[super_genes$Complex.V_v1 != ""]
  super_genes_hash[["Complex.V"]] <-
    super_genes$Complex.V_v2[super_genes$Complex.V_v2 != ""]
  
  
  state_iter = 1
  for (key_state in names(super_genes_hash)) {
    print(state_iter)
    state_iter <- state_iter + 1
    print(key_state)
    print(paste(super_genes_hash[[key_state]], sep = " "))
    
    print("+++++++++++++++")
  }
  
  
  fifteen_cluster_list_20_hash <- hash()
  co_expression_hash <- hash()
  fifteen_cluster_list_all_hash <- hash()
  new_sc_cluster_list_hash <- hash()
  new_sc_cluster_list_hash_top_hundred <- hash()
  dim(new_sc_cluster_list)[1]
  
  
  for (row in 1:nrow(fifteen_cluster_list_20)) {
    key <- as.character(fifteen_cluster_list_20[row, "cluster"])
    
    if (!has.key(key, fifteen_cluster_list_20_hash)) {
      fifteen_cluster_list_20_hash[key] <- list()
      fifteen_cluster_list_20_hash[key] <-
        append(fifteen_cluster_list_20_hash[[key]],
               fifteen_cluster_list_20[row, "gene"])
    }
    else{
      fifteen_cluster_list_20_hash[key] <-
        append(fifteen_cluster_list_20_hash[[key]],
               fifteen_cluster_list_20[row, "gene"])
    }
    
  }
  
}

# =========================================================================
# =========================================================================
# =========================================================================
# Part-2: Generating spatial seurat objects
# =========================================================================
# =========================================================================
# =========================================================================
part_2 <- TRUE
if (part_2) {
  # ----------------------------------------------- STutility 2D
  # ----------------------------------------------- STutility 2D
  # ----------------------------------------------- STutility 2D
  # ----------------------------------------------- STutility 2D
  
  samples	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-D_GTCGCGACAA-TGGCAGATTG/filtered_feature_bc_matrix.h5"
    )
  spotfiles	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-D_GTCGCGACAA-TGGCAGATTG/spatial/tissue_positions_list.csv"
    )
  imgs	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-D_GTCGCGACAA-TGGCAGATTG/spatial/tissue_hires_image.png"
    )
  json <-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-D_GTCGCGACAA-TGGCAGATTG/spatial/scalefactors_json.json"
    )
  
  infoTable <- data.frame(samples, spotfiles, imgs, json)
  se_2D <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    minSpotsPerGene = 50,
    minGenesPerSpot = 100,
    minUMICountsPerSpot = 500,
    platform =  "Visium"
  )
  
  
  head(se_2D[[]])
  p1 <-
    ST.FeaturePlot(
      se_2D,
      features = c("nFeature_RNA"),
      cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
      ncol = 2,
      pt.size = 1.3
    )
  
  str(se_2D)
  st.object <- GetStaffli(se_2D)
  st.object
  head(st.object[[]]) %>%
    kbl() %>%
    kable_styling()
  se_2D <- LoadImages(se_2D, time.resolve = FALSE, verbose = TRUE)
  ImagePlot(se_2D, method = "raster", type = "raw")
  head(se_2B[[]])
  ST.FeaturePlot(
    object = se_2D,
    features = c("nCount_RNA"),
    cols = c("lightgray", "mistyrose", "red", "dark red", "black"),
    ncol = 2,
    pt.size = 2
  )
  
  se_2D <- SCTransform(se_2D, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  se_2D <- FindNeighbors(se_2D, reduction = "pca", dims = 1:30)
  se_2D <- FindClusters(se_2D, resolution = 1.0, verbose = FALSE)
  se_2D <- RunUMAP(se_2D, reduction = "pca", dims = 1:30)
  se_2D <- SetIdent(se_2D, value = "seurat_clusters")
  FeatureOverlay(se_2D, features = "seurat_clusters",  pt.size = 2)
  # -------------------------------------------------------------------
  # ------------------------------------------------------------------- exclude stromal and immune clusters
  anchors <-
    FindTransferAnchors(
      reference = ccrcc_reference,
      query = se_2D,
      normalization.method = "SCT"
    )
  ccrcc_reference[[]]
  predictions.assay <-
    TransferData(
      anchorset = anchors,
      refdata = ccrcc_reference$annotation2,
      prediction.assay = TRUE,
      weight.reduction = se_2D[["pca"]],
      dims = 1:30
    )
  se_2D[["predictions"]] <- predictions.assay
  DefaultAssay(se_2D) <- "predictions"
  
  CellsMeta <- t(se_2D@assays$predictions$data)
  se_2D <- AddMetaData(se_2D, CellsMeta)
  
  se_2D_df <-
    se_2D[[]][c("seurat_clusters",
                "nCount_RNA",
                "nFeature_RNA",
                rownames(predictions.assay))]
  table(se_2D_df$seurat_clusters)
  se_2D_df_summary <-
    se_2D_df %>% group_by(seurat_clusters) %>% summarise_each(funs(mean))
  View(se_2D_df_summary)
  FeatureOverlay(se_2D, features = "seurat_clusters",  pt.size = 2.0)
  
  excluded_clsuters <-
    se_2D_df_summary[se_2D_df_summary[, "Urothelium-Ureter"] > 0.4 |
                       se_2D_df_summary[, "B-Cell"] > 0.5 |
                       se_2D_df_summary[, "Mesangial"] > 0.7, ]$seurat_clusters
  included_clusters <-
    setdiff(unique(se_2D[[]]$seurat_clusters), excluded_clsuters)
  
  # exclude low count clusters
  included_clusters <- setdiff(included_clusters, c("1", "5", "8"))
  #3,4,6,7,9
  tumor_2D_subset <- subset(se_2D, idents = included_clusters)
  
  FeatureOverlay(tumor_2D_subset, features = "seurat_clusters",  pt.size = 2.0)
  
  # ----------------------------------------------- STutility 2C
  # ----------------------------------------------- STutility 2C
  # ----------------------------------------------- STutility 2C
  # ----------------------------------------------- STutility 2C
  
  samples	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-C_GTTAAACCGC-CAAGCCGGCT/filtered_feature_bc_matrix.h5"
    )
  spotfiles	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-C_GTTAAACCGC-CAAGCCGGCT/spatial/tissue_positions_list.csv"
    )
  imgs	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-C_GTTAAACCGC-CAAGCCGGCT/spatial/tissue_hires_image.png"
    )
  json <-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-C_GTTAAACCGC-CAAGCCGGCT/spatial/scalefactors_json.json"
    )
  
  infoTable <- data.frame(samples, spotfiles, imgs, json)
  se_2C <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    minSpotsPerGene = 50,
    minGenesPerSpot = 100,
    minUMICountsPerSpot = 500,
    platform =  "Visium"
  )
  ST.FeaturePlot(
    se_2C,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
    ncol = 2,
    pt.size = 1.3
  )
  
  str(se_2C)
  st.object <- GetStaffli(se_2C)
  st.object
  head(st.object[[]]) %>%
    kbl() %>%
    kable_styling()
  se_2C <- LoadImages(se_2C, time.resolve = FALSE, verbose = TRUE)
  ImagePlot(se_2C, method = "raster", type = "raw")
  
  ST.FeaturePlot(
    object = se_2C,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "dark red", "black"),
    ncol = 2
  )
  
  se_2C <- SCTransform(se_2C, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  se_2C <- FindNeighbors(se_2C, reduction = "pca", dims = 1:40)
  se_2C <- FindClusters(se_2C, resolution = 1.2, verbose = FALSE)
  se_2C <- RunUMAP(se_2C, reduction = "pca", dims = 1:40)
  se_2C <- SetIdent(se_2C, value = "seurat_clusters")
  FeatureOverlay(se_2C, features = "seurat_clusters",  pt.size = 1.0)
  # -------------------------------------------------------------------
  # ------------------------------------------------------------------- exclude stromal and immune clusters
  anchors <-
    FindTransferAnchors(
      reference = ccrcc_reference,
      query = se_2C,
      normalization.method = "SCT"
    )
  ccrcc_reference[[]]
  predictions.assay <-
    TransferData(
      anchorset = anchors,
      refdata = ccrcc_reference$annotation2,
      prediction.assay = TRUE,
      weight.reduction = se_2C[["pca"]],
      dims = 1:30
    )
  se_2C[["predictions"]] <- predictions.assay
  DefaultAssay(se_2C) <- "predictions"
  
  CellsMeta <- t(se_2C@assays$predictions$data)
  se_2C <- AddMetaData(se_2C, CellsMeta)
  
  se_2C_df <-
    se_2C[[]][c("seurat_clusters",
                "nCount_RNA",
                "nFeature_RNA",
                rownames(predictions.assay))]
  table(se_2C_df$seurat_clusters)
  se_2C_df_summary <-
    se_2C_df %>% group_by(seurat_clusters) %>% summarise_each(funs(mean))
  View(se_2C_df_summary)
  FeatureOverlay(se_2C, features = "seurat_clusters",  pt.size = 2.0)
  
  excluded_clsuters <-
    se_2C_df_summary[se_2C_df_summary[, "Urothelium-Ureter"] > 0.4 |
                       se_2C_df_summary[, "B-Cell"] > 0.5 |
                       se_2C_df_summary[, "Mesangial"] > 0.7, ]$seurat_clusters
  included_clusters <-
    setdiff(unique(se_2C[[]]$seurat_clusters), excluded_clsuters)
  
  
  tumor_2C_subset <- subset(se_2C, idents = included_clusters)
  
  FeatureOverlay(tumor_2C_subset, features = "seurat_clusters",  pt.size = 2.0)
  
  # -------------------------------------------------------------------
  # -------------------------------------------------------------------
  
  #save(tumor_2C_subset, file = '~/Documents/visium_data/stu_analysis/tumor_seurat_object_2C_spatial_white.RData')
  
  
  # ----------------------------------------------- STutility 1D
  # ----------------------------------------------- STutility 1D
  # ----------------------------------------------- STutility 1D
  # ----------------------------------------------- STutility 1D
  
  samples	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-D_CGCATTTCTC-TGGTTGTGCG/filtered_feature_bc_matrix.h5"
    )
  spotfiles	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-D_CGCATTTCTC-TGGTTGTGCG/spatial/tissue_positions_list.csv"
    )
  imgs	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-D_CGCATTTCTC-TGGTTGTGCG/spatial/tissue_hires_image.png"
    )
  json <-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-D_CGCATTTCTC-TGGTTGTGCG/spatial/scalefactors_json.json"
    )
  
  infoTable <- data.frame(samples, spotfiles, imgs, json)
  se_1D <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    minSpotsPerGene = 50,
    minGenesPerSpot = 100,
    minUMICountsPerSpot = 500,
    platform =  "Visium"
  )
  ST.FeaturePlot(
    se_1D,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
    ncol = 2,
    pt.size = 1.3
  )
  
  str(se_1D)
  st.object <- GetStaffli(se_1D)
  se_1D <- LoadImages(se_1D, time.resolve = FALSE, verbose = TRUE)
  #ImagePlot(se_1D, method = "raster", type = "raw")
  
  ST.FeaturePlot(
    object = se_1D,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "dark red", "black"),
    ncol = 2
  )
  
  se_1D <- SCTransform(se_1D, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  se_1D <- FindNeighbors(se_1D, reduction = "pca", dims = 1:30)
  se_1D <- FindClusters(se_1D, resolution = 1.0, verbose = FALSE)
  se_1D <- RunUMAP(se_1D, reduction = "pca", dims = 1:30)
  se_1D <- SetIdent(se_1D, value = "seurat_clusters")
  FeatureOverlay(se_1D, features = "seurat_clusters",  pt.size = 1.0)
  # -------------------------------------------------------------------
  # ------------------------------------------------------------------- exclude stromal and immune clusters
  anchors <-
    FindTransferAnchors(
      reference = ccrcc_reference,
      query = se_1D,
      normalization.method = "SCT"
    )
  ccrcc_reference[[]]
  predictions.assay <-
    TransferData(
      anchorset = anchors,
      refdata = ccrcc_reference$annotation2,
      prediction.assay = TRUE,
      weight.reduction = se_1D[["pca"]],
      dims = 1:30
    )
  se_1D[["predictions"]] <- predictions.assay
  DefaultAssay(se_1D) <- "predictions"
  
  CellsMeta <- t(se_1D@assays$predictions$data)
  se_1D <- AddMetaData(se_1D, CellsMeta)
  
  se_1D_df <-
    se_1D[[]][c("seurat_clusters",
                "nCount_RNA",
                "nFeature_RNA",
                rownames(predictions.assay))]
  table(se_1D_df$seurat_clusters)
  se_1D_df_summary <-
    se_1D_df %>% group_by(seurat_clusters) %>% summarise_each(funs(mean))
  View(se_1D_df_summary)
  FeatureOverlay(se_1D, features = "seurat_clusters",  pt.size = 2.0)
  
  excluded_clsuters <-
    se_1D_df_summary[se_1D_df_summary[, "B-Cell"] > 0.5 |
                       se_1D_df_summary[, "Mesangial"] > 0.7, ]$seurat_clusters
  included_clusters <-
    setdiff(unique(se_1D[[]]$seurat_clusters), excluded_clsuters)
  
  # exclude low count clusters
  included_clusters <- setdiff(included_clusters, c("2"))
  #3,4,6,7,9
  tumor_1D_subset <- subset(se_1D, idents = included_clusters)
  
  FeatureOverlay(tumor_1D_subset, features = "seurat_clusters",  pt.size = 2.0)
  
  # -------------------------------------------------------------------
  # -------------------------------------------------------------------
  #save(tumor_1D_subset, file = '~/Documents/visium_data/stu_analysis/tumor_seurat_object_1D_spatial_white.RData')
  
  
  # ----------------------------------------------- STutility 1B
  # ----------------------------------------------- STutility 1B
  # ----------------------------------------------- STutility 1B
  # ----------------------------------------------- STutility 1B
  
  samples	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-B_CGAGCGCAGT-CAGTTGAGGT/filtered_feature_bc_matrix.h5"
    )
  spotfiles	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-B_CGAGCGCAGT-CAGTTGAGGT/spatial/tissue_positions_list.csv"
    )
  imgs	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-B_CGAGCGCAGT-CAGTTGAGGT/spatial/tissue_hires_image.png"
    )
  json <-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S1-B_CGAGCGCAGT-CAGTTGAGGT/spatial/scalefactors_json.json"
    )
  
  infoTable <- data.frame(samples, spotfiles, imgs, json)
  se_1B <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    minSpotsPerGene = 50,
    minGenesPerSpot = 100,
    minUMICountsPerSpot = 500,
    platform =  "Visium"
  )
  
  se_1B <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    min.gene.spots = 5,
    min.spot.count = 500,
    platform =  "Visium"
  )
  ST.FeaturePlot(
    se_1B,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
    ncol = 2,
    pt.size = 1.3
  )
  
  str(se_1B)
  st.object <- GetStaffli(se_1B)
  st.object
  head(st.object[[]]) %>%
    kbl() %>%
    kable_styling()
  se_1B <- LoadImages(se_1B, time.resolve = FALSE, verbose = TRUE)
  ImagePlot(se_1B, method = "raster", type = "raw")
  
  ST.FeaturePlot(
    object = se_1B,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "dark red", "black"),
    ncol = 2
  )
  
  
  se_1B <- SCTransform(se_1B, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  se_1B <- FindNeighbors(se_1B, reduction = "pca", dims = 1:30)
  se_1B <- FindClusters(se_1B, resolution = 1.0, verbose = FALSE)
  se_1B <- RunUMAP(se_1B, reduction = "pca", dims = 1:30)
  se_1B <- SetIdent(se_1B, value = "seurat_clusters")
  FeatureOverlay(se_1B, features = "seurat_clusters",  pt.size = 1.0)
  # -------------------------------------------------------------------
  # ------------------------------------------------------------------- exclude stromal and immune clusters
  anchors <-
    FindTransferAnchors(
      reference = ccrcc_reference,
      query = se_1B,
      normalization.method = "SCT"
    )
  ccrcc_reference[[]]
  predictions.assay <-
    TransferData(
      anchorset = anchors,
      refdata = ccrcc_reference$annotation2,
      prediction.assay = TRUE,
      weight.reduction = se_1B[["pca"]],
      dims = 1:30
    )
  se_1B[["predictions"]] <- predictions.assay
  DefaultAssay(se_1B) <- "predictions"
  
  CellsMeta <- t(se_1B@assays$predictions$data)
  se_1B <- AddMetaData(se_1B, CellsMeta)
  
  se_1B_df <-
    se_1B[[]][c("seurat_clusters",
                "nCount_RNA",
                "nFeature_RNA",
                rownames(predictions.assay))]
  table(se_1B_df$seurat_clusters)
  se_1B_df_summary <-
    se_1B_df %>% group_by(seurat_clusters) %>% summarise_each(funs(mean))
  View(se_1B_df_summary)
  FeatureOverlay(se_1B, features = "seurat_clusters",  pt.size = 2.0)
  
  excluded_clsuters <-
    se_1B_df_summary[se_1B_df_summary[, "B-Cell"] > 0.5 |
                       se_1B_df_summary[, "Mesangial"] > 0.6, ]$seurat_clusters
  included_clusters <-
    setdiff(unique(se_1B[[]]$seurat_clusters), excluded_clsuters)
  
  # exclude low count clusters
  included_clusters <- setdiff(included_clusters, c("0"))
  #3,4,6,7,9
  tumor_1B_subset <- subset(se_1B, idents = included_clusters)
  
  FeatureOverlay(tumor_1B_subset, features = "seurat_clusters",  pt.size = 2.0)
  
  # -------------------------------------------------------------------
  # -------------------------------------------------------------------
  #save(tumor_1B_subset, file = '~/Documents/visium_data/stu_analysis/tumor_seurat_object_1B_spatial_white.RData')
  
  
  # ----------------------------------------------- STutility 2B
  # ----------------------------------------------- STutility 2B
  # ----------------------------------------------- STutility 2B
  # ----------------------------------------------- STutility 2B
  # ----------------------------------------------- STutility 2B
  
  samples	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-B_GCGCAAACAA-AGTTATCTAG/filtered_feature_bc_matrix.h5"
    )
  spotfiles	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-B_GCGCAAACAA-AGTTATCTAG/spatial/tissue_positions_list.csv"
    )
  imgs	<-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-B_GCGCAAACAA-AGTTATCTAG/spatial/tissue_hires_image.png"
    )
  
  #imgs	<- c("~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-B_GCGCAAACAA-AGTTATCTAG/spatial/tissue_hires_image_white.png")
  json <-
    c(
      "~/Documents/visium_data/10x_analysis_5687-AJ/Sample_5687-AJ-S2-B_GCGCAAACAA-AGTTATCTAG/spatial/scalefactors_json.json"
    )
  
  infoTable <- data.frame(samples, spotfiles, imgs, json)
  se_2B <- InputFromTable(
    infotable = infoTable,
    minUMICountsPerGene = 100,
    minSpotsPerGene = 50,
    minGenesPerSpot = 100,
    minUMICountsPerSpot = 500,
    platform =  "Visium"
  )
  ST.FeaturePlot(
    se_2B,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "darkred", "black"),
    ncol = 2,
    pt.size = 1.3
  )
  
  str(se_2B)
  st.object <- GetStaffli(se_2B)
  st.object
  head(st.object[[]]) %>%
    kbl() %>%
    kable_styling()
  se_2B <- LoadImages(se_2B, time.resolve = FALSE, verbose = TRUE)
  ImagePlot(se_2B, method = "raster", type = "raw")
  
  ST.FeaturePlot(
    object = se_2B,
    features = c("nFeature_RNA"),
    cols = c("lightgray", "mistyrose", "red", "dark red", "black"),
    ncol = 2
  )
  
  se_2B <- SCTransform(se_2B, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  se_2B <- FindNeighbors(se_2B, reduction = "pca", dims = 1:30)
  se_2B <- FindClusters(se_2B, verbose = FALSE)
  se_2B <- RunUMAP(se_2B, reduction = "pca", dims = 1:30)
  se_2B <- SetIdent(se_2B, value = "seurat_clusters")
  
  # -------------------------------------------------------------------
  # ------------------------------------------------------------------- exclude stromal and immune clusters
  anchors <-
    FindTransferAnchors(
      reference = ccrcc_reference,
      query = se_2B,
      normalization.method = "SCT"
    )
  ccrcc_reference[[]]
  predictions.assay <-
    TransferData(
      anchorset = anchors,
      refdata = ccrcc_reference$annotation2,
      prediction.assay = TRUE,
      weight.reduction = se_2B[["pca"]],
      dims = 1:30
    )
  se_2B[["predictions"]] <- predictions.assay
  DefaultAssay(se_2B) <- "predictions"
  
  CellsMeta <- t(se_2B@assays$predictions$data)
  se_2B <- AddMetaData(se_2B, CellsMeta)
  
  se_2B_df <-
    se_2B[[]][c("seurat_clusters",
                "nCount_RNA",
                "nFeature_RNA",
                rownames(predictions.assay))]
  table(se_2B_df$seurat_clusters)
  se_2B_df_summary <-
    se_2B_df %>% group_by(seurat_clusters) %>% summarise_each(funs(mean))
  View(se_2B_df_summary)
  FeatureOverlay(se_2B, features = "seurat_clusters",  pt.size = 2.0)
  
  excluded_clsuters <-
    se_2B_df_summary[se_2B_df_summary[, "B-Cell"] > 0.5 |
                       se_2B_df_summary[, "Mesangial"] > 0.6, ]$seurat_clusters
  included_clusters <-
    setdiff(unique(se_2B[[]]$seurat_clusters), excluded_clsuters)
  #3,4,6,7,9
  tumor_2B_subset <- subset(se_2B, idents = included_clusters)
  
  FeatureOverlay(tumor_2B_subset, features = "seurat_clusters",  pt.size = 2.0)
  
}
# =========================================================================
# =========================================================================
# =========================================================================
# Part3: Selecting best resolution and integrating spatial sections
# =========================================================================
# =========================================================================
# =========================================================================
part_3 <- TRUE
if (part_3) {
  tumor_1B_subset <-
    SCTransform(tumor_1B_subset, assay = "RNA", verbose = FALSE)
  tumor_1D_subset <-
    SCTransform(tumor_1D_subset, assay = "RNA", verbose = FALSE)
  tumor_2B_subset <-
    SCTransform(tumor_2B_subset, assay = "RNA", verbose = FALSE)
  tumor_2C_subset <-
    SCTransform(tumor_2C_subset, assay = "RNA", verbose = FALSE)
  tumor_2D_subset <-
    SCTransform(tumor_2D_subset, assay = "RNA", verbose = FALSE)
  
  
  patient.list <-
    c(
      tumor_1B_subset,
      tumor_1D_subset,
      tumor_2B_subset,
      tumor_2C_subset,
      tumor_2D_subset
    )
  
  features <-
    SelectIntegrationFeatures(object.list = patient.list, nfeatures = 2000)
  patient.list <-
    PrepSCTIntegration(object.list = patient.list, anchor.features = features)
  subset.anchors <-
    FindIntegrationAnchors(
      object.list = patient.list,
      normalization.method = "SCT",
      anchor.features = features
    )
  subset.combined.sct <-
    IntegrateData(anchorset = subset.anchors, normalization.method = "SCT")
  subset.combined.sct <- RunPCA(subset.combined.sct, verbose = FALSE)
  subset.combined.sct <-
    RunUMAP(subset.combined.sct,
            reduction = "pca",
            dims = 1:30)
  subset.combined.sct <-
    FindNeighbors(subset.combined.sct,
                  reduction = "pca",
                  dims = 1:30)
  
  DefaultAssay(object = subset.combined.sct) <- "integrated"
  
  subset.combined.sct <-
    FindClusters(
      subset.combined.sct,
      resolution = c(
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
        1.0,
        1.05,
        1.1,
        1.15,
        1.2
      )
    )
  
  
  cluster_sl_eval <- function(mat, clust) {
    out = as.data.frame(matrix(rep(NA, 100 * ncol(clust)), nrow = 100))
    
    for (x in 1:100) {
      i = sample(1:ncol(mat), min(1000, ncol(mat)))
      d = as.dist(1 - cor(mat[, i], method = "pearson"))
      
      for (j in 1:ncol(clust)) {
        if (length(table(clust[i, j])) == 1) {
          out[x, j] = 0
        }
        else{
          sil = silhouette(as.numeric(clust[i, j]), d)
          out[x, j] = mean(sil[, "sil_width"])
        }
      }
    }
    means = apply(out, 2, mean)
    sd = apply(out, 2, sd)
    return(list(means, sd))
  }
  ress <-
    c(
      "integrated_snn_res.0.3",
      "integrated_snn_res.0.35",
      "integrated_snn_res.0.4",
      "integrated_snn_res.0.45",
      "integrated_snn_res.0.5",
      "integrated_snn_res.0.55",
      "integrated_snn_res.0.6",
      "integrated_snn_res.0.65",
      "integrated_snn_res.0.7",
      "integrated_snn_res.0.75",
      "integrated_snn_res.0.8",
      "integrated_snn_res.0.85",
      "integrated_snn_res.0.9",
      "integrated_snn_res.0.95",
      "integrated_snn_res.1"
    )
  ress_names <- c(
    "0.3",
    "0.35",
    "0.4",
    "0.45",
    "0.5",
    "0.55",
    "0.6",
    "0.65",
    "0.7",
    "0.75",
    "0.8",
    "0.85",
    "0.9",
    "0.95",
    "1"
  )
  
  clust = subset.combined.sct@meta.data[, ress]
  
  mat = subset.combined.sct@assays$integrated@scale.data
  
  out = cluster_sl_eval(mat, clust)
  means = out[[1]]
  sd = out[[2]]
  means
  col
  resolution_df <- data.frame(resolution  = means)
  rownames(resolution_df) <- ress_names
  #colnames(mm_decoy) <- ress
  
  setEPS()                                             # Set postscript arguments
  postscript(
    '~/Documents/visium_data/stu_analysis/plots/resolution.eps',
    width = 12,
    height = 6
  )
  
  par(mar = c(5, 8 , 4 , 1))
  #barplot(height=data$value, names=data$name, col="#69b3a2")
  
  x <- barplot(
    as.matrix(t(resolution_df)),
    beside = TRUE,
    legend.text = rownames(resolution_df),
    #args.legend=list(bty="n",horiz=TRUE),
    #width = 0.65,
    #col=color_list,border="black",
    #xlim=c(0,60),xlab="Spatial Clsuter",
    #ylim=c(0,0.12),ylab=tl2,
    main = "Resolution"
  )
  
  y <- as.matrix(resolution_df)
  #abline( h = 50, col="red", lwd=3, lty=2)
  #text(x,y,labels=as.character(y))
  dev.off()
  
  # The best resolution is: 0.55
  Idents(subset.combined.sct) <- "integrated_snn_res.0.55"
  DimPlot(subset.combined.sct,
          group.by = "integrated_snn_res.0.55",
          reduction = "umap")
  
}
# =========================================================================
# =========================================================================
# =========================================================================
# Part4: Add super gene module score for the integrated assay
# =========================================================================
# =========================================================================
# =========================================================================
part_4 <- TRUE
if (part_4) {
  glutathionine_cu_signature <-
    read.csv(file = "~/Documents/visium_data/stu_analysis/Supergene_signatures.csv", header =
               TRUE)
  
  
  for (key_state in names(glutathionine_cu_signature)) {
    print(key_state)
    print(paste(glutathionine_cu_signature[[key_state]], sep = " "))
    custom_features <-
      list(unlist(glutathionine_cu_signature[[key_state]][glutathionine_cu_signature[[key_state]] != ""]))
    print(custom_features)
  }
  
  for (key_state in names(glutathionine_cu_signature)) {
    print(key_state)
    print(glutathionine_cu_signature[[key_state]])
  }
  
  DefaultAssay(subset.combined.sct) <- "SCT"
  for (key_state in names(glutathionine_cu_signature)) {
    print(key_state)
    key_state_name <- paste0("final_integrated_", key_state)
    #print(paste0("final_integrated_",glutathionine_cu_signature[[key_state]]))
    print("key_state_name : ")
    print(key_state_name)
    custom_features <-
      list(unlist(glutathionine_cu_signature[[key_state]][glutathionine_cu_signature[[key_state]] != ""]))
    print(custom_features)
    subset.combined.sct <-  AddModuleScore(
      subset.combined.sct,
      features = custom_features,
      name = key_state_name,
      ctrl = 80,
      nbin = 25
    )
  }
  # ---------------------------------------
  lb <-
    summary(subset.combined.sct[[]]$final_integrated_Cu.signature1)["1st Qu."]
  ub <-
    summary(subset.combined.sct[[]]$final_integrated_Cu.signature1)["3rd Qu."]
  DefaultAssay(subset.combined.sct) <- "SCT"
  cu_discretized <-
    discretize_array_fun(subset.combined.sct[[]]$final_integrated_Cu.signature1,
                         lb,
                         ub)
  names(cu_discretized) <- colnames(x = subset.combined.sct)
  subset.combined.sct <- AddMetaData(object = subset.combined.sct,
                                     metadata = cu_discretized,
                                     col.name = 'Cu.signature.discretised')
  pal_bl_or <- colorRampPalette(c("blue", "lightgrey", "gold2"))
  FeaturePlot(
    subset.combined.sct,
    features = c("Cu.signature.discretised"),
    cols = pal(3)
  )
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/Cu.signature.discretised.eps')
  FeaturePlot(
    subset.combined.sct,
    features = c("Cu.signature.discretised"),
    cols = pal_bl_or(3),
    pt.size = 0.5,
    ncol = 1
  ) &
    theme(
      text = element_text(),
      axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      ),
      axis.title = element_text(size = 20),
      axis.title.y.right = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      axis.line = element_line(size = 1)
    )
  
  dev.off()
  # oxphos discretized ===========================================
  # oxphos discretized ===========================================
  # oxphos discretized ===========================================
  lb <-
    summary(subset.combined.sct[[]]$final_integrated_oxPhos1)["1st Qu."]
  ub <-
    summary(subset.combined.sct[[]]$final_integrated_oxPhos1)["3rd Qu."]
  lb
  ub
  DefaultAssay(subset.combined.sct) <- "SCT"
  cu_discretized <-
    discretize_array_fun(subset.combined.sct[[]]$final_integrated_oxPhos1, lb, ub)
  names(cu_discretized) <- colnames(x = subset.combined.sct)
  subset.combined.sct <- AddMetaData(object = subset.combined.sct,
                                     metadata = cu_discretized,
                                     col.name = 'oxphos.signature.discretised')
  FeaturePlot(
    subset.combined.sct,
    features = c("oxphos.signature.discretised"),
    cols = pal(3)
  )
  setEPS()                                             # Set postscript arguments
  postscript(
    '~/Documents/visium_data/stu_analysis/plots/oxphos.signature.discretised.eps'
  )
  FeaturePlot(
    subset.combined.sct,
    features = c("oxphos.signature.discretised"),
    cols = pal_bl_or(3),
    pt.size = 0.5,
    ncol = 1
  ) &
    theme(
      text = element_text(),
      axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      ),
      axis.title = element_text(size = 20),
      axis.title.y.right = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      axis.line = element_line(size = 1)
    )
  
  dev.off()
  # gsh discretized ===========================================
  # gsh discretized ===========================================
  # gsh discretized ===========================================
  lb <-
    summary(subset.combined.sct[[]]$final_integrated_Glutathione1)["1st Qu."]
  ub <-
    summary(subset.combined.sct[[]]$final_integrated_Glutathione1)["3rd Qu."]
  lb
  ub
  DefaultAssay(subset.combined.sct) <- "SCT"
  cu_discretized <-
    discretize_array_fun(subset.combined.sct[[]]$final_integrated_Glutathione1, lb, ub)
  names(cu_discretized) <- colnames(x = subset.combined.sct)
  subset.combined.sct <- AddMetaData(object = subset.combined.sct,
                                     metadata = cu_discretized,
                                     col.name = 'gsh.signature.discretised')
  FeaturePlot(
    subset.combined.sct,
    features = c("gsh.signature.discretised"),
    cols = pal(3)
  )
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/gsh.signature.discretised.eps')
  FeaturePlot(
    subset.combined.sct,
    features = c("gsh.signature.discretised"),
    cols = pal_bl_or(3),
    pt.size = 0.5,
    ncol = 1
  ) &
    theme(
      text = element_text(),
      axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      ),
      axis.title = element_text(size = 20),
      axis.title.y.right = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      axis.line = element_line(size = 1)
    )
  
  dev.off()
  # ---
}
# =========================================================================
# =========================================================================
# =========================================================================
# Part5: Comparison of module scores of super genes
# =========================================================================
# =========================================================================
# =========================================================================
part_5 <- TRUE
if (part_5) {
  query_pathways_names <-
    c(
      "0",
      "1",
      "10",
      "11",
      "12",
      "13",
      "2",
      "4",
      "5",
      "6"
      ,
      "7",
      "8",
      "9",
      "Cu",
      "EMT"
      ,
      "FAO",
      "FAS",
      "glycolysis",
      "HIF1A",
      "HIF1A.2A"
      ,
      "MAS",
      "MTs",
      "NRF2_targets",
      "oxPhos",
      "Glutathione"
      ,
      "Cu.II",
      "Cu.III",
      "Cu.signature"
    )
  col_annot_all_bins <- query_pathways_names
  
  row_annot_all_bins <- query_pathways_names
  
  #load('~/Documents/visium_data/stu_analysis/query_pathways.RData')
  #query_pathways:
  #[1] "all_bins_01"            "all_bins_11"            "all_bins_101"           "all_bins_111"           "all_bins_121"
  #[6] "all_bins_131"           "all_bins_21"            "all_bins_41"            "all_bins_51"            "all_bins_61"
  #[11] "all_bins_71"            "all_bins_81"            "all_bins_91"            "all_bins_Cu1"           "all_bins_EMT1"
  #[16] "all_bins_FAO1"          "all_bins_FAS1"          "all_bins_glycolysis1"   "all_bins_HIF1A1"        "all_bins_HIF1A.2A1"
  #[21] "all_bins_MAS1"          "all_bins_MTs1"          "all_bins_NRF2_targets1" "all_bins_oxPhos1"       "Glutathione1"
  #[26] "Cu.II1"                 "Cu.III1"                "Cu.signature1"
  data_subset_1B_all_bins = matrix(0, 28, 28)
  data_subset_1D_all_bins = matrix(0, 28, 28)
  data_subset_2B_all_bins = matrix(0, 28, 28)
  data_subset_2C_all_bins = matrix(0, 28, 28)
  data_subset_2D_all_bins = matrix(0, 28, 28)
  
  data_subset_1B_all_bins_p_val = matrix(0, 28, 28)
  data_subset_1D_all_bins_p_val = matrix(0, 28, 28)
  data_subset_2B_all_bins_p_val = matrix(0, 28, 28)
  data_subset_2C_all_bins_p_val = matrix(0, 28, 28)
  data_subset_2D_all_bins_p_val = matrix(0, 28, 28)
  colnames(tumor_1D_subset[[]])
  # ================================== 1D
  
  
  row_iter <- 0
  for (titem in query_pathways) {
    col_iter <- 0
    row_iter <- row_iter + 1
    for (titem2 in query_pathways) {
      col_iter <- col_iter + 1
      
      v1 <- tumor_1D_subset[[]][, titem]
      v2 <- tumor_1D_subset[[]][, titem2]
      
      res <- cor.test(v1, v2,  method = "pearson")
      
      data_subset_1D_all_bins[row_iter, col_iter] <- res$estimate
      data_subset_1D_all_bins_p_val[row_iter, col_iter] <-
        res$p.value
      
      #print("----------------")
    }
  }
  row_annot_all_bins
  rownames(data_subset_1D_all_bins) = row_annot_all_bins
  colnames(data_subset_1D_all_bins) = col_annot_all_bins
  rownames(data_subset_1D_all_bins_p_val) = row_annot_all_bins
  colnames(data_subset_1D_all_bins_p_val) = col_annot_all_bins
  png(
    file = "~/Documents/visium_data/stu_analysis/heatmap/Slide_1D_all_bins_heatmap.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 300
  )
  
  Heatmap(data_subset_1D_all_bins,
          column_title = "Slide 1D_all_bins")
  dev.off()
  
  # ================================== 1B
  
  row_iter <- 0
  for (titem in query_pathways) {
    col_iter <- 0
    row_iter <- row_iter + 1
    for (titem2 in query_pathways) {
      col_iter <- col_iter + 1
      #print(tumor_1D_subset[[]][titem])
      #print(tumor_1D_subset[[]][titem2])
      v1 <- tumor_1B_subset[[]][, titem]
      v2 <- tumor_1B_subset[[]][, titem2]
      
      res <- cor.test(v1, v2,  method = "pearson")
      
      data_subset_1B_all_bins[row_iter, col_iter] <- res$estimate
      data_subset_1B_all_bins_p_val[row_iter, col_iter] <-
        res$p.value
      #print("----------------")
    }
  }
  rownames(data_subset_1B_all_bins) = row_annot_all_bins
  colnames(data_subset_1B_all_bins) = col_annot_all_bins
  rownames(data_subset_1B_all_bins_p_val) = row_annot_all_bins
  colnames(data_subset_1B_all_bins_p_val) = col_annot_all_bins
  png(
    file = "~/Documents/visium_data/stu_analysis/heatmap/Slide_1B_all_bins_heatmap.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 300
  )
  
  Heatmap(data_subset_1B_all_bins,
          column_title = "Slide 1B_all_bins")
  dev.off()
  # ================================== 2B
  row_iter <- 0
  for (titem in query_pathways) {
    col_iter <- 0
    row_iter <- row_iter + 1
    for (titem2 in query_pathways) {
      col_iter <- col_iter + 1
      #print(tumor_1D_subset[[]][titem])
      #print(tumor_1D_subset[[]][titem2])
      v1 <- tumor_2B_subset[[]][, titem]
      v2 <- tumor_2B_subset[[]][, titem2]
      
      res <- cor.test(v1, v2,  method = "pearson")
      
      data_subset_2B_all_bins[row_iter, col_iter] <- res$estimate
      data_subset_2B_all_bins_p_val[row_iter, col_iter] <-
        res$p.value
      #print("----------------")
    }
  }
  rownames(data_subset_2B_all_bins) = row_annot_all_bins
  colnames(data_subset_2B_all_bins) = col_annot_all_bins
  rownames(data_subset_2B_all_bins_p_val) = row_annot_all_bins
  colnames(data_subset_2B_all_bins_p_val) = col_annot_all_bins
  png(
    file = "~/Documents/visium_data/stu_analysis/heatmap/Slide_2B_all_bins_heatmap.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 300
  )
  
  Heatmap(data_subset_2B_all_bins,
          column_title = "Slide 2B_all_bins")
  dev.off()
  # ================================== 2C
  row_iter <- 0
  for (titem in query_pathways) {
    col_iter <- 0
    row_iter <- row_iter + 1
    for (titem2 in query_pathways) {
      col_iter <- col_iter + 1
      #print(tumor_1D_subset[[]][titem])
      #print(tumor_1D_subset[[]][titem2])
      v1 <- tumor_2C_subset[[]][, titem]
      v2 <- tumor_2C_subset[[]][, titem2]
      
      res <- cor.test(v1, v2,  method = "pearson")
      
      data_subset_2C_all_bins[row_iter, col_iter] <- res$estimate
      data_subset_2C_all_bins_p_val[row_iter, col_iter] <-
        res$p.value
      #print("----------------")
    }
  }
  rownames(data_subset_2C_all_bins) = row_annot_all_bins
  colnames(data_subset_2C_all_bins) = col_annot_all_bins
  rownames(data_subset_2C_all_bins_p_val) = row_annot_all_bins
  colnames(data_subset_2C_all_bins_p_val) = col_annot_all_bins
  png(
    file = "~/Documents/visium_data/stu_analysis/heatmap/Slide_2C_all_bins_heatmap.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 300
  )
  
  Heatmap(data_subset_2C_all_bins,
          column_title = "Slide 2C_all_bins")
  dev.off()
  # ================================== 2D
  row_iter <- 0
  for (titem in query_pathways) {
    col_iter <- 0
    row_iter <- row_iter + 1
    for (titem2 in query_pathways) {
      col_iter <- col_iter + 1
      #print(tumor_1D_subset[[]][titem])
      #print(tumor_1D_subset[[]][titem2])
      v1 <- tumor_2D_subset[[]][, titem]
      v2 <- tumor_2D_subset[[]][, titem2]
      
      res <- cor.test(v1, v2,  method = "pearson")
      
      data_subset_2D_all_bins[row_iter, col_iter] <- res$estimate
      data_subset_2D_all_bins_p_val[row_iter, col_iter] <-
        res$p.value
      #print("----------------")
    }
  }
  rownames(data_subset_2D_all_bins) = row_annot_all_bins
  colnames(data_subset_2D_all_bins) = col_annot_all_bins
  rownames(data_subset_2D_all_bins_p_val) = row_annot_all_bins
  colnames(data_subset_2D_all_bins_p_val) = col_annot_all_bins
  png(
    file = "~/Documents/visium_data/stu_analysis/heatmap/Slide_2D_all_bins_heatmap.png",
    width = 10,
    height = 10,
    units = 'in',
    res = 300
  )
  
  Heatmap(data_subset_2D_all_bins,
          column_title = "Slide 2D_all_bins")
  dev.off()
  
  # --------------------------------------------
  # glut/cu heatmap --------------------------
  # --------------------------------------------
  heatmap_df_pval <-
    data.frame(matrix(0,    # Create empty data frame
                      nrow = 5,
                      ncol = 7))
  
  rownames(heatmap_df_pval) <- c("2C", "2D", "1B", "1D", "2B")
  
  colnames(heatmap_df_pval) <-
    c(
      'DF/RL',
      'Glutathione/Cu.signature',
      'Glutathione/OxPhos',
      'Glutathione/NRF2',
      'Cu.signature/OxPhos',
      'Cu.signature/NRF2',
      'OxPhos/NRF2'
    )
  heatmap_df_pval$'DF/RL' <- c("DF", "DF", "RL", "RL", "RL")
  
  heatmap_df_pval["2B", 'Glutathione/Cu.signature'] <-
    data_subset_2B_all_bins_p_val["Glutathione", "Cu.signature"]
  heatmap_df_pval["2B", 'Glutathione/OxPhos'] <-
    data_subset_2B_all_bins_p_val["Glutathione", "oxPhos"]
  heatmap_df_pval["2B", 'Glutathione/NRF2'] <-
    data_subset_2B_all_bins_p_val["Glutathione", "NRF2_targets"]
  heatmap_df_pval["2B", 'Cu.signature/OxPhos'] <-
    data_subset_2B_all_bins_p_val["oxPhos", "Cu.signature"]
  heatmap_df_pval["2B", 'Cu.signature/NRF2'] <-
    data_subset_2B_all_bins_p_val["NRF2_targets", "Cu.signature"]
  heatmap_df_pval["2B", 'OxPhos/NRF2'] <-
    data_subset_2B_all_bins_p_val["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df_pval["2C", 'Glutathione/Cu.signature'] <-
    data_subset_2C_all_bins_p_val["Glutathione", "Cu.signature"]
  heatmap_df_pval["2C", 'Glutathione/OxPhos'] <-
    data_subset_2C_all_bins_p_val["Glutathione", "oxPhos"]
  heatmap_df_pval["2C", 'Glutathione/NRF2'] <-
    data_subset_2C_all_bins_p_val["Glutathione", "NRF2_targets"]
  heatmap_df_pval["2C", 'Cu.signature/OxPhos'] <-
    data_subset_2C_all_bins_p_val["oxPhos", "Cu.signature"]
  heatmap_df_pval["2C", 'Cu.signature/NRF2'] <-
    data_subset_2C_all_bins_p_val["NRF2_targets", "Cu.signature"]
  heatmap_df_pval["2C", 'OxPhos/NRF2'] <-
    data_subset_2C_all_bins_p_val["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df_pval["2D", 'Glutathione/Cu.signature'] <-
    data_subset_2D_all_bins_p_val["Glutathione", "Cu.signature"]
  heatmap_df_pval["2D", 'Glutathione/OxPhos'] <-
    data_subset_2D_all_bins_p_val["Glutathione", "oxPhos"]
  heatmap_df_pval["2D", 'Glutathione/NRF2'] <-
    data_subset_2D_all_bins_p_val["Glutathione", "NRF2_targets"]
  heatmap_df_pval["2D", 'Cu.signature/OxPhos'] <-
    data_subset_2D_all_bins_p_val["oxPhos", "Cu.signature"]
  heatmap_df_pval["2D", 'Cu.signature/NRF2'] <-
    data_subset_2D_all_bins_p_val["NRF2_targets", "Cu.signature"]
  heatmap_df_pval["2D", 'OxPhos/NRF2'] <-
    data_subset_2D_all_bins_p_val["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df_pval["1B", 'Glutathione/Cu.signature'] <-
    data_subset_1B_all_bins_p_val["Glutathione", "Cu.signature"]
  heatmap_df_pval["1B", 'Glutathione/OxPhos'] <-
    data_subset_1B_all_bins_p_val["Glutathione", "oxPhos"]
  heatmap_df_pval["1B", 'Glutathione/NRF2'] <-
    data_subset_1B_all_bins_p_val["Glutathione", "NRF2_targets"]
  heatmap_df_pval["1B", 'Cu.signature/OxPhos'] <-
    data_subset_1B_all_bins_p_val["Cu.signature", "oxPhos"]
  heatmap_df_pval["1B", 'Cu.signature/NRF2'] <-
    data_subset_1B_all_bins_p_val["NRF2_targets", "Cu.signature"]
  heatmap_df_pval["1B", 'OxPhos/NRF2'] <-
    data_subset_1B_all_bins_p_val["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df_pval["1D", 'Glutathione/Cu.signature'] <-
    data_subset_1D_all_bins_p_val["Glutathione", "Cu.signature"]
  heatmap_df_pval["1D", 'Glutathione/OxPhos'] <-
    data_subset_1D_all_bins_p_val["Glutathione", "oxPhos"]
  heatmap_df_pval["1D", 'Glutathione/NRF2'] <-
    data_subset_1D_all_bins_p_val["Glutathione", "NRF2_targets"]
  heatmap_df_pval["1D", 'Cu.signature/OxPhos'] <-
    data_subset_1D_all_bins_p_val["oxPhos", "Cu.signature"]
  heatmap_df_pval["1D", 'Cu.signature/NRF2'] <-
    data_subset_1D_all_bins_p_val["NRF2_targets", "Cu.signature"]
  heatmap_df_pval["1D", 'OxPhos/NRF2'] <-
    data_subset_1D_all_bins_p_val["NRF2_targets", "oxPhos"]
  
  #=============================================
  heatmap_df <- data.frame(matrix(0,    # Create empty data frame
                                  nrow = 5,
                                  ncol = 7))
  
  rownames(heatmap_df) <- c("2C", "2D", "1B", "1D", "2B")
  
  colnames(heatmap_df) <-
    c(
      'DF/RL',
      'Glutathione/Cu.signature',
      'Glutathione/OxPhos',
      'Glutathione/NRF2',
      'Cu.signature/OxPhos',
      'Cu.signature/NRF2',
      'OxPhos/NRF2'
    )
  heatmap_df$'DF/RL' <- c("DF", "DF", "RL", "RL", "RL")
  colnames(data_subset_2B_all_bins)
  heatmap_df["2B", 'Glutathione/Cu.signature'] <-
    data_subset_2B_all_bins["Glutathione", "Cu.signature"]
  heatmap_df["2B", 'Glutathione/OxPhos'] <-
    data_subset_2B_all_bins["Glutathione", "oxPhos"]
  heatmap_df["2B", 'Glutathione/NRF2'] <-
    data_subset_2B_all_bins["Glutathione", "NRF2_targets"]
  heatmap_df["2B", 'Cu.signature/OxPhos'] <-
    data_subset_2B_all_bins["oxPhos", "Cu.signature"]
  heatmap_df["2B", 'Cu.signature/NRF2'] <-
    data_subset_2B_all_bins["NRF2_targets", "Cu.signature"]
  heatmap_df["2B", 'OxPhos/NRF2'] <-
    data_subset_2B_all_bins["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df["2C", 'Glutathione/Cu.signature'] <-
    data_subset_2C_all_bins["Glutathione", "Cu.signature"]
  heatmap_df["2C", 'Glutathione/OxPhos'] <-
    data_subset_2C_all_bins["Glutathione", "oxPhos"]
  heatmap_df["2C", 'Glutathione/NRF2'] <-
    data_subset_2C_all_bins["Glutathione", "NRF2_targets"]
  heatmap_df["2C", 'Cu.signature/OxPhos'] <-
    data_subset_2C_all_bins["oxPhos", "Cu.signature"]
  heatmap_df["2C", 'Cu.signature/NRF2'] <-
    data_subset_2C_all_bins["NRF2_targets", "Cu.signature"]
  heatmap_df["2C", 'OxPhos/NRF2'] <-
    data_subset_2C_all_bins["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df["2D", 'Glutathione/Cu.signature'] <-
    data_subset_2D_all_bins["Glutathione", "Cu.signature"]
  heatmap_df["2D", 'Glutathione/OxPhos'] <-
    data_subset_2D_all_bins["Glutathione", "oxPhos"]
  heatmap_df["2D", 'Glutathione/NRF2'] <-
    data_subset_2D_all_bins["Glutathione", "NRF2_targets"]
  heatmap_df["2D", 'Cu.signature/OxPhos'] <-
    data_subset_2D_all_bins["oxPhos", "Cu.signature"]
  heatmap_df["2D", 'Cu.signature/NRF2'] <-
    data_subset_2D_all_bins["NRF2_targets", "Cu.signature"]
  heatmap_df["2D", 'OxPhos/NRF2'] <-
    data_subset_2D_all_bins["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df["1B", 'Glutathione/Cu.signature'] <-
    data_subset_1B_all_bins["Glutathione", "Cu.signature"]
  heatmap_df["1B", 'Glutathione/OxPhos'] <-
    data_subset_1B_all_bins["Glutathione", "oxPhos"]
  heatmap_df["1B", 'Glutathione/NRF2'] <-
    data_subset_1B_all_bins["Glutathione", "NRF2_targets"]
  heatmap_df["1B", 'Cu.signature/OxPhos'] <-
    data_subset_1B_all_bins["Cu.signature", "oxPhos"]
  heatmap_df["1B", 'Cu.signature/NRF2'] <-
    data_subset_1B_all_bins["NRF2_targets", "Cu.signature"]
  heatmap_df["1B", 'OxPhos/NRF2'] <-
    data_subset_1B_all_bins["NRF2_targets", "oxPhos"]
  
  #------------------
  
  heatmap_df["1D", 'Glutathione/Cu.signature'] <-
    data_subset_1D_all_bins["Glutathione", "Cu.signature"]
  heatmap_df["1D", 'Glutathione/OxPhos'] <-
    data_subset_1D_all_bins["Glutathione", "oxPhos"]
  heatmap_df["1D", 'Glutathione/NRF2'] <-
    data_subset_1D_all_bins["Glutathione", "NRF2_targets"]
  heatmap_df["1D", 'Cu.signature/OxPhos'] <-
    data_subset_1D_all_bins["oxPhos", "Cu.signature"]
  heatmap_df["1D", 'Cu.signature/NRF2'] <-
    data_subset_1D_all_bins["NRF2_targets", "Cu.signature"]
  heatmap_df["1D", 'OxPhos/NRF2'] <-
    data_subset_1D_all_bins["NRF2_targets", "oxPhos"]
  
  #------------------
  #------------------
  #------------------
  Glutathione_Cu.signature <- heatmap_df$`Glutathione.Cu.signature`
  Glutathione_OxPhos <- heatmap_df$`Glutathione.OxPhos`
  #Glutathione_NRF2 <- heatmap_df$`Glutathione.NRF2`
  Cu.signature_OxPhos <- heatmap_df$`Cu.signature.OxPhos`
  #Cu.signature_NRF2 <- heatmap_df$`Cu.signature.NRF2`
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/glut_copper_oxphos_final.eps')
  #eps('~/Documents/visium_data/stu_analysis/metabolic_Cu_glu_Non_overlapping_glu_box_plot.tiff', width = 5, height = 3, units = 'in', res = 300)
  par(cex.main = 2)
  
  par(cex.lab = 2) # is for y-axis
  par(cex.axis = 2) # is for x-axis
  
  #OxPhos_NRF2
  c1 <-
    match.closest(mean(Glutathione_Cu.signature), pal_21_num, tolerance = 0.1)
  c2 <-
    match.closest(mean(Cu.signature_OxPhos), pal_21_num, tolerance = 0.1)
  
  c3 <-
    match.closest(mean(Glutathione_OxPhos), pal_21_num, tolerance = 0.1)
  
  
  boxplot(
    Glutathione_Cu.signature,
    Cu.signature_OxPhos,
    Glutathione_OxPhos,
    
    
    main = "module comparision",
    at = c(1, 2, 3),
    names = c("Cu/GSH", "OxPhos/Cu", "OxPhos/GSH"),
    las = 2,
    boxwex = 0.5,
    col = c(pal_21[c1], pal_21[c2], pal_21[c3]),
    
    border = "brown",
    horizontal = TRUE,
    notch = FALSE,
    pars = list(par(mar = c(7, 14, 4, 2)))
  )
  title(xlab = "Correlation Coefficient",
        cex.lab = 2,
        line = 4.5)
  text(
    x = median(Glutathione_Cu.signature),
    y = 1.4,
    labels = paste("p <", format(
      mean(Glutathione_Cu.signature_pval), scientific = TRUE
    )),
    cex = 1.2
  )
  text(
    x = median(Cu.signature_OxPhos),
    y = 2.4,
    labels = paste("p <", format(
      mean(Cu.signature_OxPhos_pval), scientific = TRUE
    )),
    cex = 1.2
  )
  text(
    x = median(Glutathione_OxPhos),
    y = 3.4,
    labels = paste("p <", format(
      mean(Glutathione_OxPhos_pval), scientific = TRUE
    )),
    cex = 1.2
  )
  #text(x= median(Glutathione_NRF2), y= 7.6, labels= paste("p-value <", format(max(Glutathione_NRF2_pval), scientific = TRUE)))
  dev.off()
  
  
}
# =========================================================================
# =========================================================================
# =========================================================================
# Part6: Analysis of single cell to spatial cluster similarity
# =========================================================================
# =========================================================================
# =========================================================================
part_6 <- TRUE
if (part_6) {
  sc_clusters <-
    c(
      "z_scored_Cluster_01",
      "z_scored_Cluster_11",
      "z_scored_Cluster_101",
      "z_scored_Cluster_111"
      ,
      "z_scored_Cluster_121",
      "z_scored_Cluster_131",
      "z_scored_Cluster_141",
      "z_scored_Cluster_21"
      ,
      "z_scored_Cluster_31",
      "z_scored_Cluster_41",
      "z_scored_Cluster_51",
      "z_scored_Cluster_61"
      ,
      "z_scored_Cluster_71",
      "z_scored_Cluster_81",
      "z_scored_Cluster_91"
    )
  
  sc_clusters <-
    c(
      "Cluster_01",
      "Cluster_11",
      "Cluster_101",
      "Cluster_111",
      "Cluster_121",
      "Cluster_131",
      "Cluster_141",
      "Cluster_21",
      "Cluster_31",
      "Cluster_41",
      "Cluster_51",
      "Cluster_61",
      "Cluster_71",
      "Cluster_81",
      "Cluster_91"
    )
  
  sp_clusters <- c("sp_0",
                   "sp_1",
                   "sp_2",
                   "sp_3",
                   "sp_4",
                   "sp_5",
                   "sp_6",
                   "sp_7",
                   "sp_8")
  sc_clusters_names <-
    c(
      "sc_0",
      "sc_1",
      "sc_10",
      "sc_11",
      "sc_12",
      "sc_13",
      "sc_14",
      "sc_2",
      "sc_3",
      "sc_4",
      "sc_5",
      "sc_6",
      "sc_7",
      "sc_8",
      "sc_9"
    )
  
  all_distributions <- c()
  sp_sc_distribution <- hash()
  
  for (sp in sp_clusters) {
    sp_sc_distribution[sp] <- hash()
    for (sc in sc_clusters_names) {
      sp_sc_distribution[[sp]][sc] <- c()
    }
  }
  
  colnames(tumor_1B_subset[[]])
  for (obj in c(
    tumor_1B_subset,
    tumor_1D_subset,
    tumor_2B_subset,
    tumor_2C_subset,
    tumor_2D_subset
  )) {
    #print(obj)
    print("-----")
    for (col in 1:length(sc_clusters)) {
      col_name <- sc_clusters[col]
      print(col_name)
      scaled_col <- scale(obj[[]][, col_name])
      
      iter <- 1
      #print(obj[[]][,col_name])
      for (row in rownames(obj[[]])) {
        #print(row)
        row_name <- paste0("sp_", obj[[]][row, "seurat_clusters"])
        
        sp_sc_distribution[[row_name]][[sc_clusters_names[col]]] <-
          append(sp_sc_distribution[[row_name]][[sc_clusters_names[col]]], scaled_col[iter])
        
        
        iter <- iter + 1
        
      }
      
    }
  }
  
  
  second_corr2 <- matrix(0.0, 9, 15)
  rownames(second_corr2) = sp_clusters
  colnames(second_corr2) = sc_clusters_names
  
  row_corr_df = c()
  col_corr_df = c()
  p_val_corr_df = c()
  mean_corr_d = c()
  for (sp in sp_clusters) {
    sp_new_name <- gsub("_", ":", sp)
    
    for (sc in sc_clusters_names) {
      sc_new_name <- gsub("_", ":", sc)
      print(sc)
      print(sp)
      
      row_corr_df <- append(row_corr_df, sp_new_name)
      col_corr_df <- append(col_corr_df, sc_new_name)
      mean_corr_d <-
        append(mean_corr_d, mean(sp_sc_distribution[[sp]][[sc]]))
      
      tt <- t.test(sp_sc_distribution[[sp]][[sc]], mu = 0)
      
      
      if (tt$p.value == 0) {
        p_val_corr_df <- append(p_val_corr_df, 200)
      } else{
        p_val_corr_df <- append(p_val_corr_df, -log10(tt$p.value))
        
      }
      
      
    }
    print("--------")
  }
  
  corr_df <-
    data.frame(
      sp = row_corr_df,
      sc = col_corr_df,
      p_val = p_val_corr_df,
      mean = mean_corr_d
    )
  
  cut_slots <-
    seq(max(corr_df$mean) + 0.2, min(corr_df$mean) - 0.2,  by = -0.2)
  corr_df$mean_cut <- cut(corr_df$mean, breaks = cut_slots)
  summary(corr_df$mean_cut)
  options(repr.plot.width = 8, repr.plot.height = 4)
  
  pal_yel_11_mod <- pal_yel(11)
  
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/corr_new.eps')
  
  ggplot(corr_df,
         aes(x = col_corr_df, y = row_corr_df),
         vjust = 0.3,
         hjust = 0.5) +
    geom_tile(color = "black", fill = "white") +
    geom_point(aes(color = mean_cut, size = abs(p_val))) +
    
    scale_color_manual(values = pal_yel_11_mod) +
    scale_size_continuous(range = c(2, 10), name = "-log10(p)") + theme(aspect.ratio =
                                                                          0.6)
  
  dev.off()
  
  
}

# =========================================================================
# =========================================================================
# =========================================================================
# Part7: Cluster localization analysis
# =========================================================================
# =========================================================================
# =========================================================================
part_7 <- TRUE
if (part_7) {
  ne_an_1B <-
    neighborhood_histogram_analysis(tumor_1B_subset, neighbor_df_1B)
  dim(ne_an_1B)
  decoy_clusters <- ne_an_1B$class_decoy
  names(decoy_clusters) <- colnames(x = tumor_1B_subset)
  tumor_1B_subset <- AddMetaData(object = tumor_1B_subset,
                                 metadata = decoy_clusters,
                                 col.name = 'decoy_clusters')
  colnames(tumor_2D_subset[[]])
  # ==============================
  ne_an_1D <-
    neighborhood_histogram_analysis(tumor_1D_subset, neighbor_df_1D)
  decoy_clusters <- ne_an_1D$class_decoy
  names(decoy_clusters) <- colnames(x = tumor_1D_subset)
  tumor_1D_subset <- AddMetaData(object = tumor_1D_subset,
                                 metadata = decoy_clusters,
                                 col.name = 'decoy_clusters')
  # ==============================
  ne_an_2B <-
    neighborhood_histogram_analysis(tumor_2B_subset, neighbor_df_2B)
  decoy_clusters <- ne_an_2B$class_decoy
  names(decoy_clusters) <- colnames(x = tumor_2B_subset)
  tumor_2B_subset <- AddMetaData(object = tumor_2B_subset,
                                 metadata = decoy_clusters,
                                 col.name = 'decoy_clusters')
  # ==============================
  ne_an_2D <-
    neighborhood_histogram_analysis(tumor_2D_subset, neighbor_df_2D)
  decoy_clusters <- ne_an_2D$class_decoy
  names(decoy_clusters) <- colnames(x = tumor_2D_subset)
  tumor_2D_subset <- AddMetaData(object = tumor_2D_subset,
                                 metadata = decoy_clusters,
                                 col.name = 'decoy_clusters')
  # ==============================
  ne_an_2C <-
    neighborhood_histogram_analysis(tumor_2C_subset, neighbor_df_2C)
  decoy_clusters <- ne_an_2C$class_decoy
  names(decoy_clusters) <- colnames(x = tumor_2C_subset)
  tumor_2C_subset <- AddMetaData(object = tumor_2C_subset,
                                 metadata = decoy_clusters,
                                 col.name = 'decoy_clusters')
  # ==============================
  ne_an <- rbind(ne_an_1B, ne_an_1D)
  ne_an <- rbind(ne_an, ne_an_2B)
  ne_an <- rbind(ne_an, ne_an_2C)
  ne_an <- rbind(ne_an,ne_an_2D)
  
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/neighbor_histo_nov.eps')
  ggplot(ne_an, aes(x = value, y = class, fill = stat(x))) + geom_density_ridges_gradient() + scale_fill_viridis_c(name = "Depth", option = "C")
  dev.off()
  
  
  head(ne_an)
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/neighbor_histo_decoy_nov.eps')
  ggplot(ne_an, aes(x = value_decoy, y = class_decoy, fill = stat(x))) + geom_density_ridges_gradient() + scale_fill_viridis_c(name = "Depth", option = "C")
  dev.off()
  
  setEPS()                                             # Set postscript arguments
  postscript('~/Documents/visium_data/stu_analysis/plots/neighbor_histo_1d.eps')
  ggplot(ne_an_1D, aes(x = value, y = class, fill = stat(x))) + geom_density_ridges_gradient()
  + scale_fill_viridis_c(name = "Depth", option = "C")#+ theme(axis.text = element_text(size = 20),text = element_text(size = 20))
  dev.off()
}


# =========================================================================
# =========================================================================
# =========================================================================
# Part8: Cluster neighborhood analysis
# =========================================================================
# =========================================================================
# =========================================================================
part_8 <- TRUE
if (part_8) {
  subset.combined.sct_1B <-
    subset(subset.combined.sct, subset = slide == '1B')
  subset.combined.sct_1D <-
    subset(subset.combined.sct, subset = slide == '1D')
  subset.combined.sct_2B <-
    subset(subset.combined.sct, subset = slide == '2B')
  subset.combined.sct_2C <-
    subset(subset.combined.sct, subset = slide == '2C')
  subset.combined.sct_2D <-
    subset(subset.combined.sct, subset = slide == '2D')
  colnames(tumor_2B_subset[[]])
  #==========================================
  # 1B
  #==========================================
  ccrcc_sobject_small <- RenameCells(object = subset.combined.sct_1B,
                            new.names = substr(colnames(x = subset.combined.sct_1B), 1, nchar(colnames(x = subset.combined.sct_1B)) -
                                                 2))
  
  CellsMeta = ccrcc_sobject_small@meta.data
  CellsMetaTrim <- subset(CellsMeta, select = c("seurat_clusters"))
  #CellsMetaTrim <- subset(CellsMeta, select = c("integrated_snn_res.0.8"))
  colnames(CellsMetaTrim) <- "seurat_clusters"
  tumor_1B_subset <- AddMetaData(tumor_1B_subset, CellsMetaTrim)
  summary(subset.combined.sct_1B$seurat_clusters, 20)
  head(tumor_1B_subset$seurat_clusters, 20)
  #==========================================
  # 1D
  #==========================================
  ccrcc_sobject_small <- RenameCells(object = subset.combined.sct_1D,
                            new.names = substr(colnames(x = subset.combined.sct_1D), 1, nchar(colnames(x = subset.combined.sct_1D)) -
                                                 2))
  
  CellsMeta = ccrcc_sobject_small@meta.data
  CellsMetaTrim <- subset(CellsMeta, select = c("seurat_clusters"))
  #CellsMetaTrim <- subset(CellsMeta, select = c("integrated_snn_res.0.8"))
  colnames(CellsMetaTrim) <- "seurat_clusters"
  tumor_1D_subset <- AddMetaData(tumor_1D_subset, CellsMetaTrim)
  head(subset.combined.sct_1D$seurat_clusters, 20)
  head(tumor_1D_subset$seurat_clusters, 20)
  #==========================================
  # 2B
  #==========================================
  ccrcc_sobject_small <- RenameCells(object = subset.combined.sct_2B,
                            new.names = substr(colnames(x = subset.combined.sct_2B), 1, nchar(colnames(x = subset.combined.sct_2B)) -
                                                 2))
  
  CellsMeta = ccrcc_sobject_small@meta.data
  CellsMetaTrim <- subset(CellsMeta, select = c("seurat_clusters"))
  #CellsMetaTrim <- subset(CellsMeta, select = c("integrated_snn_res.0.8"))
  colnames(CellsMetaTrim) <- "seurat_clusters"
  tumor_2B_subset <- AddMetaData(tumor_2B_subset, CellsMetaTrim)
  head(subset.combined.sct_2B$seurat_clusters, 20)
  head(tumor_2B_subset$seurat_clusters, 20)
  #==========================================
  # 2C
  #==========================================
  ccrcc_sobject_small <- RenameCells(object = subset.combined.sct_2C,
                            new.names = substr(colnames(x = subset.combined.sct_2C), 1, nchar(colnames(x = subset.combined.sct_2C)) -
                                                 2))
  
  CellsMeta = ccrcc_sobject_small@meta.data
  CellsMetaTrim <- subset(CellsMeta, select = c("seurat_clusters"))
  #CellsMetaTrim <- subset(CellsMeta, select = c("integrated_snn_res.0.8"))
  colnames(CellsMetaTrim) <- "seurat_clusters"
  tumor_2C_subset <- AddMetaData(tumor_2C_subset, CellsMetaTrim)
  head(subset.combined.sct_2C$seurat_clusters, 20)
  head(tumor_2C_subset$seurat_clusters, 20)
  #==========================================
  # 2D
  #==========================================
  ccrcc_sobject_small <- RenameCells(object = subset.combined.sct_2D,
                            new.names = substr(colnames(x = subset.combined.sct_2D), 1, nchar(colnames(x = subset.combined.sct_2D)) -
                                                 2))
  
  CellsMeta = ccrcc_sobject_small@meta.data
  CellsMetaTrim <- subset(CellsMeta, select = c("seurat_clusters"))
  #CellsMetaTrim <- subset(CellsMeta, select = c("integrated_snn_res.0.8"))
  colnames(CellsMetaTrim) <- "seurat_clusters"
  head(CellsMetaTrim)
  
  tumor_2D_subset <- AddMetaData(tumor_2D_subset, CellsMetaTrim)
  head(subset.combined.sct_2D$seurat_clusters, 20)
  head(tumor_2D_subset$seurat_clusters, 20)
  colnames(subset.combined.sct_2D[[]])
  # ======================================
  plot_zero_vs_others_neighbor <- TRUE
  if (plot_zero_vs_others_neighbor) {
    ne_an_1B_for_one_cluster <-
      neighborhood_histogram_analysis_for_one_cluster_reverse(tumor_1B_subset, neighbor_df_1B, 0)
    dim(ne_an_1B_for_one_cluster)
    table(ne_an_1B_for_one_cluster$class_decoy)
    decoy_clusters <- as.numeric(ne_an_1B_for_one_cluster$class_decoy)
    names(decoy_clusters) <- colnames(x = tumor_1B_subset)
    tumor_1B_subset <- AddMetaData(object = tumor_1B_subset,
                                   metadata = decoy_clusters,
                                   col.name = 'decoy_clusters_fixed_zero')
    tumor_1B_subset[[]]$decoy_clusters_fixed_zero
    tumor_1B_subset[[]]$seurat_clusters
    head(ne_an_1B_for_one_cluster)
    head(tumor_1B_subset[[]]$decoy_clusters_fixed_zero)
    head(tumor_1B_subset[[]]$seurat_clusters)
    FeatureOverlay(
      tumor_1B_subset,
      features = c("decoy_clusters_fixed_zero"),
      cols = color_list
    )
    FeatureOverlay(tumor_1B_subset,
                   features = c("seurat_clusters"),
                   cols = color_list)
    
    
    ne_an_1D_for_one_cluster <-
      neighborhood_histogram_analysis_for_one_cluster_reverse(tumor_1D_subset_final, neighbor_df_1D, 0)
    ne_an_2B_for_one_cluster <-
      neighborhood_histogram_analysis_for_one_cluster_reverse(tumor_2B_subset_final, neighbor_df_2B, 0)
    ne_an_2C_for_one_cluster <-
      neighborhood_histogram_analysis_for_one_cluster_reverse(tumor_2C_subset, neighbor_df_2C, 0)
    ne_an_2D_for_one_cluster <-
      neighborhood_histogram_analysis_for_one_cluster_reverse(tumor_2D_subset_final, neighbor_df_2D, 0)
    "seurat_clusters" %in% colnames(tumor_1D_subset_final[[]])
    total_classes_counts_list <-
      tumor_1D_subset_final[[]][, "seurat_clusters"]
    total_classes_counts_list <-
      append(total_classes_counts_list, tumor_1B_subset[[]][, "seurat_clusters"])
    total_classes_counts_list <-
      append(total_classes_counts_list, tumor_2B_subset_final[[]][, "seurat_clusters"])
    total_classes_counts_list <-
      append(total_classes_counts_list, tumor_2C_subset[[]][, "seurat_clusters"])
    tumor_1D_subset_final[[]][, "seurat_clusters"]
    total_classes_counts <- table(total_classes_counts_list)
    total_classes_counts[8]
    ne_an_2 <- rbind(ne_an_1B_for_one_cluster, ne_an_1D_for_one_cluster)
    
    ne_an_2 <- rbind(ne_an_2, ne_an_2B_for_one_cluster)
    
    ne_an_2 <- rbind(ne_an_2, ne_an_2C_for_one_cluster)
    
    # ------------------------------
    # ------------------------------
    # ------------------------------
    
    
    tail(ne_an_2)
    ne_an_2_vectorized1 <-
      data.frame(class = rep(, desired_length),
                 value = ne_an_2$orig_1)
    # ------------------------------
    # ------------------------------
    # ------------------------------
    # ------------------------------
    stack_plot_neigh <- c()
    stack_plot_number <- c()
    stack_plot_cluster <- c()
    stack_plot_class <- c()
    # ------------------------------
    # ------------------------------
    # ------------------------------111
    ne_an_1_orig <- ne_an_2$orig_1
    
    P <- table(ne_an_1_orig)
    
    ne_an_1_decoy <- ne_an_2$decoy_1
    table(ne_an_1_decoy)
    Q <- table(ne_an_1_decoy)
    
    P["0"] <-
      total_classes_counts[2] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[2] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    
    #Q["6"] <- 0
    #Q["5"] <- 0
    P
    PP <- P / total_classes_counts[2]
    QQ <- Q / total_classes_counts[2]
    library(philentropy)
    P
    Q
    
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-1")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-1")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP1 <- matrix(0.0, 2, 2)
    spatial_corrP1[1, 1] <- P[1] + P[2]
    spatial_corrP1[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP1[1, 2] <- Q[1] + Q[2]
    spatial_corrP1[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP1 <- matrix(0.0, 2, 2)
    spatial_corrPP1[1, 1] <- PP[1] + PP[2]
    spatial_corrPP1[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP1[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP1[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    if (FALSE) {
      spatial_corrP1 <- matrix(0.0, 2, 2)
      spatial_corrP1[1, 1] <- P[1]
      spatial_corrP1[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP1[1, 2] <- Q[1]
      spatial_corrP1[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP1 <- matrix(0.0, 2, 2)
      spatial_corrPP1[1, 1] <- PP[1]
      spatial_corrPP1[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP1[1, 2] <- QQ[1]
      spatial_corrPP1[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP1
    chisq.test(spatial_corrP1, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("1"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("1"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    length(P)
    value_list <- value_list_iter
    class_list <- class_list_iter
    value_decoy_list <- value_decoy_list_iter
    class_decoy_list <- class_decoy_list_iter
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_1_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[2] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 1" ,
      xlab = "number of neighbors in cluster 1",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[2] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    
    
    
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------222
    ne_an_2_orig <- ne_an_2$orig_2
    P <- table(ne_an_2_orig)
    
    ne_an_2_decoy <- ne_an_2$decoy_2
    table(ne_an_2_decoy)
    Q <- table(ne_an_2_decoy)
    P["0"] <-
      total_classes_counts[3] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[3] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    P
    Q
    
    PP <- P / total_classes_counts[3]
    
    QQ <- Q / total_classes_counts[3]
    library(philentropy)
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-2")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-2")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP2 <- matrix(0.0, 2, 2)
    spatial_corrP2[1, 1] <- P[1] + P[2]
    spatial_corrP2[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP2[1, 2] <- Q[1] + Q[2]
    spatial_corrP2[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP2 <- matrix(0.0, 2, 2)
    spatial_corrPP2[1, 1] <- PP[1] + PP[2]
    spatial_corrPP2[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP2[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP2[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    spatial_corrP2
    if (FALSE) {
      spatial_corrP2 <- matrix(0.0, 2, 2)
      spatial_corrP2[1, 1] <- P[1]
      spatial_corrP2[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP2[1, 2] <- Q[1]
      spatial_corrP2[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP2 <- matrix(0.0, 2, 2)
      spatial_corrPP2[1, 1] <- PP[1]
      spatial_corrPP2[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP2[1, 2] <- QQ[1]
      spatial_corrPP2[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP2
    chisq.test(spatial_corrP2, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("2"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("2"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_2_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[3] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 2" ,
      xlab = "number of neighbors in cluster 2",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[3] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------333
    ne_an_3_orig <- as.data.frame(lapply(ne_an_2["orig_3"], unlist))
    P <- table(ne_an_3_orig)
    
    ne_an_3_decoy <- as.data.frame(lapply(ne_an_2["decoy_3"], unlist))
    table(ne_an_3_decoy)
    Q <- table(ne_an_3_decoy)
    P
    Q
    
    P["6"] <- 0
    P["0"] <-
      total_classes_counts[4] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[4] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    PP <- P / total_classes_counts[4]
    QQ <- Q / total_classes_counts[4]
    library(philentropy)
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-3")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-3")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP3 <- matrix(0.0, 2, 2)
    spatial_corrP3[1, 1] <- P[1] + P[2]
    spatial_corrP3[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP3[1, 2] <- Q[1] + Q[2]
    spatial_corrP3[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP3 <- matrix(0.0, 2, 2)
    spatial_corrPP3[1, 1] <- PP[1] + PP[2]
    spatial_corrPP3[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP3[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP3[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    spatial_corrP3
    if (FALSE) {
      spatial_corrP3 <- matrix(0.0, 2, 2)
      spatial_corrP3[1, 1] <- P[1]
      spatial_corrP3[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP3[1, 2] <- Q[1]
      spatial_corrP3[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP3 <- matrix(0.0, 2, 2)
      spatial_corrPP3[1, 1] <- PP[1]
      spatial_corrPP3[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP3[1, 2] <- QQ[1]
      spatial_corrPP3[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP3
    chisq.test(spatial_corrP3, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("3"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("3"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_3_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[4] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 3",
      xlab = "number of neighbors in cluster 3",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[4] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------444
    ne_an_4_orig <- as.data.frame(lapply(ne_an_2["orig_4"], unlist))
    P <- table(ne_an_4_orig)
    
    ne_an_4_decoy <- as.data.frame(lapply(ne_an_2["decoy_4"], unlist))
    table(ne_an_4_decoy)
    Q <- table(ne_an_4_decoy)
    P["0"] <-
      total_classes_counts[5] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[5] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    P
    Q
    PP <- P / total_classes_counts[5]
    QQ <- Q / total_classes_counts[5]
    library(philentropy)
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-4")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-4")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP4 <- matrix(0.0, 2, 2)
    spatial_corrP4[1, 1] <- P[1] + P[2]
    spatial_corrP4[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP4[1, 2] <- Q[1] + Q[2]
    spatial_corrP4[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP4 <- matrix(0.0, 2, 2)
    spatial_corrPP4[1, 1] <- PP[1] + PP[2]
    spatial_corrPP4[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP4[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP4[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    
    if (FALSE) {
      spatial_corrP4 <- matrix(0.0, 2, 2)
      spatial_corrP4[1, 1] <- P[1]
      spatial_corrP4[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP4[1, 2] <- Q[1]
      spatial_corrP4[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP4 <- matrix(0.0, 2, 2)
      spatial_corrPP4[1, 1] <- PP[1]
      spatial_corrPP4[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP4[1, 2] <- QQ[1]
      spatial_corrPP4[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP4
    chisq.test(spatial_corrP4, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("4"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("4"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_4_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[5] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 4" ,
      xlab = "number of neighbors in cluster 4",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[5] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------555
    ne_an_5_orig <- as.data.frame(lapply(ne_an_2["orig_5"], unlist))
    P <- table(ne_an_5_orig)
    
    ne_an_5_decoy <- as.data.frame(lapply(ne_an_2["decoy_5"], unlist))
    table(ne_an_5_decoy)
    Q <- table(ne_an_5_decoy)
    
    P
    Q
    
    P["5"] <- 0
    P["6"] <- 0
    P["0"] <-
      total_classes_counts[6] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[6] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    PP <- P / total_classes_counts[6]
    QQ <- Q / total_classes_counts[6]
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-5")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-5")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP5 <- matrix(0.0, 2, 2)
    spatial_corrP5[1, 1] <- P[1] + P[2]
    spatial_corrP5[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP5[1, 2] <- Q[1] + Q[2]
    spatial_corrP5[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP5 <- matrix(0.0, 2, 2)
    spatial_corrPP5[1, 1] <- PP[1] + PP[2]
    spatial_corrPP5[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP5[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP5[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    if (FALSE) {
      spatial_corrP5 <- matrix(0.0, 2, 2)
      spatial_corrP5[1, 1] <- P[1]
      spatial_corrP5[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP5[1, 2] <- Q[1]
      spatial_corrP5[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP5 <- matrix(0.0, 2, 2)
      spatial_corrPP5[1, 1] <- PP[1]
      spatial_corrPP5[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP5[1, 2] <- QQ[1]
      spatial_corrPP5[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP5
    chisq.test(spatial_corrP5, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("5"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("5"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_5_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[6] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 5" ,
      xlab = "number of neighbors in cluster 5",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[6] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------666
    ne_an_6_orig <- as.data.frame(lapply(ne_an_2["orig_6"], unlist))
    P <- table(ne_an_6_orig)
    
    ne_an_6_decoy <- as.data.frame(lapply(ne_an_2["decoy_6"], unlist))
    table(ne_an_6_decoy)
    Q <- table(ne_an_6_decoy)
    library(philentropy)
    P
    Q
    
    
    Q["6"] <- 0
    P["0"] <-
      total_classes_counts[7] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[7] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    PP <- P / total_classes_counts[7]
    QQ <- Q / total_classes_counts[7]
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-6")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-6")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP6 <- matrix(0.0, 2, 2)
    spatial_corrP6[1, 1] <- P[1] + P[2]
    spatial_corrP6[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP6[1, 2] <- Q[1] + Q[2]
    spatial_corrP6[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP6 <- matrix(0.0, 2, 2)
    spatial_corrPP6[1, 1] <- PP[1] + PP[2]
    spatial_corrPP6[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP6[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP6[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    if (FALSE) {
      spatial_corrP6 <- matrix(0.0, 2, 2)
      spatial_corrP6[1, 1] <- P[1]
      spatial_corrP6[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP6[1, 2] <- Q[1]
      spatial_corrP6[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP6 <- matrix(0.0, 2, 2)
      spatial_corrPP6[1, 1] <- PP[1]
      spatial_corrPP6[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP6[1, 2] <- QQ[1]
      spatial_corrPP6[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP6
    chisq.test(spatial_corrP6, correct = FALSE)
    
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("6"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("6"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_6_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[7], rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 6",
      xlab = "number of neighbors in cluster 6",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[7] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------777
    ne_an_7_orig <- as.data.frame(lapply(ne_an_2["orig_7"], unlist))
    P <- table(ne_an_7_orig)
    
    ne_an_7_decoy <- as.data.frame(lapply(ne_an_2["decoy_7"], unlist))
    table(ne_an_7_decoy)
    Q <- table(ne_an_7_decoy)
    P
    Q
    Q["6"] <- 0
    P["0"] <-
      total_classes_counts[8] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[8] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    PP <- P / total_classes_counts[8]
    QQ <- Q / total_classes_counts[8]
    
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-7")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-7")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP7 <- matrix(0.0, 2, 2)
    spatial_corrP7[1, 1] <- P[1] + P[2]
    spatial_corrP7[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP7[1, 2] <- Q[1] + Q[2]
    spatial_corrP7[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP7 <- matrix(0.0, 2, 2)
    spatial_corrPP7[1, 1] <- PP[1] + PP[2]
    spatial_corrPP7[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP7[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP7[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    if (FALSE) {
      spatial_corrP7 <- matrix(0.0, 2, 2)
      spatial_corrP7[1, 1] <- P[1]
      spatial_corrP7[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP7[1, 2] <- Q[1]
      spatial_corrP7[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP7 <- matrix(0.0, 2, 2)
      spatial_corrPP7[1, 1] <- PP[1]
      spatial_corrPP7[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP7[1, 2] <- QQ[1]
      spatial_corrPP7[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP7
    chisq.test(spatial_corrP7, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("7"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("7"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_7_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[8] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 7",
      xlab = "number of neighbors in cluster 7",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[8] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    
    # ------------------------------
    # ------------------------------
    # ------------------------------888
    ne_an_8_orig <- as.data.frame(lapply(ne_an_2["orig_8"], unlist))
    
    P <- table(ne_an_8_orig)
    
    ne_an_8_decoy <- as.data.frame(lapply(ne_an_2["decoy_8"], unlist))
    table(ne_an_8_decoy)
    Q <- table(ne_an_8_decoy)
    P
    Q
    
    
    Q["6"] <- 0
    P["0"] <-
      total_classes_counts[9] - P[2] - P[3] - P[4] - P[5] - P[6] - P[7]
    Q["0"] <-
      total_classes_counts[9] - Q[2] - Q[3] - Q[4] - Q[5] - Q[6] - Q[7]
    PP <- P / total_classes_counts[9]
    QQ <- Q / total_classes_counts[9]
    P
    Q
    PP
    QQ
    
    for (i in 1:7)
    {
      stack_plot_neigh <- append(stack_plot_neigh, i - 1)
      
      stack_plot_number <- append(stack_plot_number, PP[i] - QQ[i])
      if (PP[i] - QQ[i] > 0) {
        stack_plot_cluster <- append(stack_plot_cluster, "sp-8")
      } else{
        stack_plot_cluster <- append(stack_plot_cluster, "negative")
      }
      stack_plot_class <- append(stack_plot_class, "sp-8")
    }
    # =========================
    # =========================
    # =========================
    # =========================
    x <- rbind(PP, QQ)
    KL(x, unit = 'log')
    ks.test(PP, QQ)
    # =========================
    # =========================
    # =========================
    # =========================
    spatial_corrP8 <- matrix(0.0, 2, 2)
    spatial_corrP8[1, 1] <- P[1] + P[2]
    spatial_corrP8[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7]
    spatial_corrP8[1, 2] <- Q[1] + Q[2]
    spatial_corrP8[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7]
    
    spatial_corrPP8 <- matrix(0.0, 2, 2)
    spatial_corrPP8[1, 1] <- PP[1] + PP[2]
    spatial_corrPP8[2, 1] <- PP[3] + PP[4] + PP[5] + PP[6] + PP[7]
    spatial_corrPP8[1, 2] <- QQ[1] + QQ[2]
    spatial_corrPP8[2, 2] <- QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7]
    if (FALSE) {
      spatial_corrP8 <- matrix(0.0, 2, 2)
      spatial_corrP8[1, 1] <- P[1]
      spatial_corrP8[2, 1] <- P[3] + P[4] + P[5] + P[6] + P[7] + P[2]
      spatial_corrP8[1, 2] <- Q[1]
      spatial_corrP8[2, 2] <- Q[3] + Q[4] + Q[5] + Q[6] + Q[7] + Q[2]
      
      spatial_corrPP8 <- matrix(0.0, 2, 2)
      spatial_corrPP8[1, 1] <- PP[1]
      spatial_corrPP8[2, 1] <-
        PP[3] + PP[4] + PP[5] + PP[6] + PP[7] + PP[2]
      spatial_corrPP8[1, 2] <- QQ[1]
      spatial_corrPP8[2, 2] <-
        QQ[3] + QQ[4] + QQ[5] + QQ[6] + QQ[7] + QQ[2]
    }
    spatial_corrP8
    chisq.test(spatial_corrP8, correct = FALSE)
    # =====================================================================================
    # =====================================================================================
    value_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(P["0"], P["1"], P["2"], P["3"], P["4"], P["5"], P["6"]))
    class_list_iter <-
      rep(c("8"), times = c(P["0"] + P["1"] + P["2"] + P["3"] + P["4"] + P["5"] +
                              P["6"]))
    value_decoy_list_iter <-
      rep(c(0, 1, 2, 3, 4, 5, 6), times = c(Q["0"], Q["1"], Q["2"], Q["3"], Q["4"], Q["5"], Q["6"]))
    class_decoy_list_iter <-
      rep(c("8"), times = c(Q["0"] + Q["1"] + Q["2"] + Q["3"] + Q["4"] + Q["5"] +
                              Q["6"]))
    
    value_list <- append(value_list, value_list_iter)
    class_list <- append(class_list, class_list_iter)
    value_decoy_list <-
      append(value_decoy_list, value_decoy_list_iter)
    class_decoy_list <-
      append(class_decoy_list, class_decoy_list_iter)
    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    data <- data.frame(
      name = c("0", "0", "1", "1", "2", "2" , "3", "3", "4", "4" , "5", "5", "6", "6"),
      average = c(PP[1], QQ[1], PP[2], QQ[2], PP[3], QQ[3], PP[4], QQ[4], PP[5], QQ[5], PP[6], QQ[6], PP[7], QQ[7]),
      number = c(P[1], Q[1], P[2], Q[2], P[3], Q[3], P[4], Q[4], P[5], Q[5], P[6], Q[6], P[7], Q[7])
    )
    P
    P[1]
    Q
    # Increase bottom margin
    par(mar = c(6, 4, 4, 4))
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/0_8_only_zero_reversed.eps',
             width = 7,
             height = 3)
    
    my_bar <- barplot(
      data$average ,
      border = F ,
      names.arg = data$name ,
      las = 2 ,
      col = c(color_list[9] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      ylim = c(0, 1) ,
      main = "Neighborhood analysis of cluster 0 with cluster 8",
      xlab = "number of neighbors in cluster 8",
      ylab = "percentage"
    )
    
    #text(my_bar, data$average+0.05 , paste("n:", data$number, sep="") ,cex=1)
    
    legend(
      "topright",
      legend = c("Slide Spots", "Shuffled Spots") ,
      col = c(color_list[9] , rgb(0.3, 0.5, 0.4, 0.6)) ,
      bty = "n",
      pch = 20 ,
      pt.cex = 2,
      cex = 0.8,
      horiz = FALSE,
      inset = c(0.05, 0.05)
    )
    
    dev.off()
    # ----------------------------------------
    # ----------------------------------------
    # ----------------------------------------
    # ----------------------------------------
    chisq.test(spatial_corrP1, correct = FALSE)
    chisq.test(spatial_corrP2, correct = FALSE)
    chisq.test(spatial_corrP3, correct = FALSE)
    chisq.test(spatial_corrP4, correct = FALSE)
    # ----------------------------------------
    chisq.test(spatial_corrP5, correct = FALSE)
    chisq.test(spatial_corrP6, correct = FALSE)
    chisq.test(spatial_corrP7, correct = FALSE)
    chisq.test(spatial_corrP8, correct = FALSE)
    # ----------------------------------------
    
    spatial_corr <- matrix(0.0, 2, 16)
    
    # ---------------------------------------
    spatial_corr[1, 1] <- spatial_corrPP7[1, 1]
    spatial_corr[2, 1] <- spatial_corrPP7[2, 1]
    spatial_corr[1, 2] <- spatial_corrPP7[1, 2]
    spatial_corr[2, 2] <- spatial_corrPP7[2, 2]
    # ---------------------------------------
    spatial_corr[1, 3] <- spatial_corrPP3[1, 1]
    spatial_corr[2, 3] <- spatial_corrPP3[2, 1]
    spatial_corr[1, 4] <- spatial_corrPP3[1, 2]
    spatial_corr[2, 4] <- spatial_corrPP3[2, 2]
    # ---------------------------------------
    spatial_corr[1, 5] <- spatial_corrPP5[1, 1]
    spatial_corr[2, 5] <- spatial_corrPP5[2, 1]
    spatial_corr[1, 6] <- spatial_corrPP5[1, 2]
    spatial_corr[2, 6] <- spatial_corrPP5[2, 2]
    # ---------------------------------------
    spatial_corr[1, 7] <- spatial_corrPP1[1, 1]
    spatial_corr[2, 7] <- spatial_corrPP1[2, 1]
    spatial_corr[1, 8] <- spatial_corrPP1[1, 2]
    spatial_corr[2, 8] <- spatial_corrPP1[2, 2]
    # ---------------------------------------
    spatial_corr[1, 9] <- spatial_corrPP4[1, 1]
    spatial_corr[2, 9] <- spatial_corrPP4[2, 1]
    spatial_corr[1, 10] <- spatial_corrPP4[1, 2]
    spatial_corr[2, 10] <- spatial_corrPP4[2, 2]
    # ---------------------------------------
    spatial_corr[1, 11] <- spatial_corrPP8[1, 1]
    spatial_corr[2, 11] <- spatial_corrPP8[2, 1]
    spatial_corr[1, 12] <- spatial_corrPP8[1, 2]
    spatial_corr[2, 12] <- spatial_corrPP8[2, 2]
    # ---------------------------------------
    spatial_corr[1, 13] <- spatial_corrPP6[1, 1]
    spatial_corr[2, 13] <- spatial_corrPP6[2, 1]
    spatial_corr[1, 14] <- spatial_corrPP6[1, 2]
    spatial_corr[2, 14] <- spatial_corrPP6[2, 2]
    # ---------------------------------------
    spatial_corr[1, 15] <- spatial_corrPP2[1, 1]
    spatial_corr[2, 15] <- spatial_corrPP2[2, 1]
    spatial_corr[1, 16] <- spatial_corrPP2[1, 2]
    spatial_corr[2, 16] <- spatial_corrPP2[2, 2]
    spatial_corr
    pal_yel <- colorRampPalette(c("blue", "gold2", "grey"))
    pal_yel_21 <- pal_yel(5)
    pal_yel <- colorRampPalette(c("darkred", "red"))
    pal_yel_21 <- pal_yel(5)
    
    sp_clusters_new <- c("no neighbors", "one or more neighbors")
    sc_clusters_names_new <-
      c(
        "sp:7",
        "sp:7_shuffled",
        "sp:3",
        "sp:3_shuffled",
        "sp:5",
        "sp:5_shuffled"
        ,
        "sp:1",
        "sp:1_shuffled",
        "sp:4",
        "sp:4_shuffled",
        "sp:6",
        "sp:6_shuffled"
        ,
        "sp:8",
        "sp:8_shuffled",
        "sp:2",
        "sp:2_shuffled"
      )
    rownames(spatial_corr) = sp_clusters_new
    colnames(spatial_corr) = sc_clusters_names_new
    setEPS()                                             # Set postscript arguments
    postscript(
      '~/Documents/visium_data/stu_analysis/plots/spatial_corr_neighbor_for_0_blue.eps'
    )
    
    corrplot(
      spatial_corr_0_1neighbor,
      method = "circle",
      is.corr = FALSE,
      col = pal_yel_21
    )
    dev.off()
    # ----------------------------------------
    # ----------------------------------------
    # ----------------------------------------
    # ----------------------------------------
    spatial_corr_0_2neighbor <- spatial_corr
    
    spatial_corr_0_1neighbor <- spatial_corr
    
    library(ggplot2)
    library(dplyr)
    library(scales)
    #color_list
    color_list_stacked <- color_list
    color_list_stacked[1] <- "darkgrey"
    cols_stack <-
      c(
        "sp-1" = color_list[2],
        "sp-2" = color_list[3],
        "sp-3" = color_list[4],
        "sp-4" = color_list[5],
        "sp-5" = color_list[6],
        "sp-6" = color_list[7],
        "sp-7" = color_list[8],
        "sp-8" = color_list[8],
        "negative" = "darkgrey"
      )
    
    stack_plot_df <-
      data.frame(
        number = stack_plot_number,
        cluster = stack_plot_cluster,
        nneighbor = stack_plot_neigh,
        class = stack_plot_class
      )
    
    cairo_ps(filename = '~/Documents/visium_data/stu_analysis/plots/stacked_barplots3.eps',
             width = 7,
             height = 3)
    stack_plot_df |>
      ggplot(aes(
        x = nneighbor,
        y = 100 * number ,
        fill = class
      ))  + scale_fill_manual(values = color_list_stacked) +
      
      geom_col() +
      facet_grid( ~ class) +
      coord_flip() #+
    #ylim(c(-25, 25))
    #geom_text(aes(label = round(number, 1)), hjust = -.1)
    dev.off()
    
    
  }
  
  
  
}





