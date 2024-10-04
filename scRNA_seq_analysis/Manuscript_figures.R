# Load libraries
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(scclusteval)
library(stringr)
library(monocle3)
library(RColorBrewer)
library(infercnv)
library(future)
setwd("./scRNA_seq_analysis/")
source("Pipelines.R")
source("Infercnv.R")
source("tools.R")
source("Hyperparameters_for_pipeline.R")
scRNA_obj = readRDS(
  paste0(intermediate_dir,"columbia_british_individually_harmony_ccRCC_harmony.rds")
)

# define target resolutions
target_res = 0.53
target_res_para = paste0("RNA_snn_res.", target_res)
Idents(scRNA_obj) = scRNA_obj@meta.data[,target_res_para]
# define signatures from file and composite expression calculation tools
source("Initialize_signature.R")
source("Signature_expression_tools.R")

stage_dic = c("S1", "S2", "S3", "Me")
names(stage_dic) = c("Stage1", "Stage2", "Stage3", "Metastatic")
scRNA_obj@meta.data[, "stage"] = stage_dic[scRNA_obj@meta.data[, "stage"]]
# find all available stages in the cohort
library(stringr)
looping_stages = as.character(sort(factor(
  unique(scRNA_obj@meta.data[, "stage"]),
  levels = stage_dic
)))

#initialize signature list
signature_list = list(
  "HIF" = signature_list_master$HIF1A,
  "Glycolysis" = signature_list_master$Glycolysis,
  "OxPhos" = c(
    signature_list_master$Complex_I,
    signature_list_master$Complex_II,
    signature_list_master$Complex_III,
    signature_list_master$Complex_IV,
    signature_list_master$Complex_V
  ),
  "Cu" = signature_list_master$`Cu-II_short`,
  "GSH" = signature_list_master$Glutathione,
  "NRF2" = signature_list_master$NRF2
)
# convert gene name from lower case to upper
signature_list = lapply(signature_list, toupper)
# extract names for this signature list
signature_list_names = names(signature_list)

# to profile the expression pattern of pathway signatures,
# we use the original RNA assay
DefaultAssay(scRNA_obj) = "RNA"
scRNA_obj = NormalizeData(scRNA_obj)

# whether remove MT high clusters or use batch corrected data to profile the expression
use_batch_corrected_data = F
remove_cluster5and6 = F
assay_in_use = "RNA"

if (remove_cluster5and6) {
  # remove cluster 5 and 6 due to their high MT expression
  cell.use = rownames(scRNA_obj@meta.data[(!scRNA_obj@meta.data[, target_res_para] %in% c(5, 6)), ])
  scRNA_obj = scRNA_obj[, cell.use]
  if (use_batch_corrected_data) {
    assay_in_use = "integrated_seurat"
    # import batch corrected data
    scRNA_obj_batch_corrected = readRDS("../scRNA_data/columbia_british_individually_harmony_ccRCC_Seurat_2000genes_plus_4supergenes.rds")
    scRNA_obj_batch_corrected = scRNA_obj_batch_corrected[, cell.use]
    # add an assay for batch corrected data
    scRNA_obj[["RNA"]]@data = as(scRNA_obj_batch_corrected[["integrated"]]@data, "dgCMatrix")
  }
} else{
  if (use_batch_corrected_data) {
    assay_in_use = "integrated_seurat"
    # inject batch corrected data
    scRNA_obj_batch_corrected = readRDS("../scRNA_data/columbia_british_individually_harmony_ccRCC_Seurat_2000genes_plus_4supergenes.rds")
    batch_corrected_data <- GetAssayData(object =  scRNA_obj_batch_corrected[['integrated']], slot = 'data')
    scRNA_obj[["integrated_seurat"]] <- CreateAssayObject(data = batch_corrected_data)
    DefaultAssay(scRNA_obj) = "integrated_seurat"
  }
}

# Add module score for for all pathway signatures
scRNA_obj = AddModuleScore(scRNA_obj,
                           features = signature_list,
                           assay = assay_in_use,
                           nbin = 24,
                           ctrl = 100,
                           name = signature_list_names)
scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)] =
  scRNA_obj@meta.data[, paste0(signature_list_names, 1:length(signature_list_names))]

# Add z-score transformed module score
if(length(signature_list_names)==1){
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")] = 
    scale(scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)])
}else{
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")] =
    as.data.frame(lapply(scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)], scale))
}



#remove redundant features
remove_column_index = match(paste0(signature_list_names, 1:length(signature_list_names)),
                            names(scRNA_obj@meta.data))
scRNA_obj@meta.data =
  scRNA_obj@meta.data[, -remove_column_index]


break_cutoffs = c(-Inf,-1.5,-0.5, 0.5, 1.5, Inf)
label_tags = create_label_tag(break_cutoffs)
if (length(signature_list_names) == 1) {
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_group")] =
    sapply(
      scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")],
      FUN = function(vx) {
        vx = cut(vx, breaks = break_cutoffs, labels = label_tags)
      }
    )
} else{
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_group")] =
    lapply(
      scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")],
      FUN = function(vx) {
        vx = cut(vx, breaks = break_cutoffs, labels = label_tags)
      }
    )
}

#######################################################
# Umaps
#######################################################
signature_list_umap = list(
  "Glycolysis" = signature_list_master$Glycolysis,
  "OxPhos" = c(
    signature_list_master$Complex_I,
    signature_list_master$Complex_II,
    signature_list_master$Complex_III,
    signature_list_master$Complex_IV,
    signature_list_master$Complex_V
  ),
  "Cu" = signature_list_master$`Cu-II_short`,
  "GSH" = signature_list_master$Glutathione
)
# convert gene name from lower case to upper
signature_list_umap = lapply(signature_list_umap, toupper)
# extract names for this signature list
signature_list_umap_names = names(signature_list_umap)

plot_title_size = 15
axis_title_size = 15
legend_text_size = 12
legend_shape_size = 5
axis_line_size = 0.5
axis_text_size = 12
umap_label_size = 1.5



# loop each stage and get umaps individually
target_signature_list_names = signature_list_umap_names
subpopulations_umpas = list()
feature_plots = list()
#define the dummy plot var to avoid changing of original data
plot_var = "target"
pmt = 0
pmr = 0
pmb = 0
pml = 0

ggplotColours <- function(n = 6, h = c(0, 360) + 15){  
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n  
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

data_bar_plot_new <- ggplotColours(n=15)

data_bar_plot_new[2] <- "darkolivegreen"
data_bar_plot_new[3] <- "cornflowerblue"
data_bar_plot_new[4] <- "cyan"
data_bar_plot_new[5] <- "gold"
data_bar_plot_new[7] <- "darkorange"
data_bar_plot_new[12] <- "goldenrod"
data_bar_plot_new[11] <- "red"
data_bar_plot_new[8] <- "blue"
data_bar_plot_new[9] <- "purple"
data_bar_plot_new[10] <- "green"


###############################################################################
###############################################################################
###############################################################################
# Figure 7A, 7B and 7C
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Convert Seurat object to Monocle 3 CellDataSet
ccRCC_obj <- scRNA_obj
count_data <- ccRCC_obj@assays$RNA@counts
cell_metadata <- ccRCC_obj@meta.data
gene_metadata <- ccRCC_obj@assays$RNA@meta.features
cds <- new_cell_data_set(
  expression_data = count_data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Define target resolution
target_res = 0.53
target_res_para = paste0("RNA_snn_res.", target_res)

# Preprocessing of trajectory analysis
recreate.partition = rep(1, length(cds@colData@rownames))
names(recreate.partition) = cds@colData@rownames
recreate.partition = as.factor(recreate.partition)
cds@clusters$UMAP$partitions = recreate.partition
Idents(ccRCC_obj) = ccRCC_obj@meta.data[,target_res_para]
cds@clusters$UMAP$clusters = ccRCC_obj@active.ident
cds@int_colData@listData$reducedDims$UMAP = ccRCC_obj@reductions$umap@cell.embeddings

# Learn the trajectory graph
cds = learn_graph(cds, use_partition = F)

# Order pseudo time (Note: When you run this line of code, there will be a window 
# that pops out to let you choose the starting nodes of the trajectory. If you want
# to reproduce the results in the paper, please select nodes from the top right where
# subclone 11 locates and then click "Done". After click "Done", you can then run
# the next line of code)
cds <- order_cells(cds, reduction_method = "UMAP")

# Visualize results
p1 = DimPlot(ccRCC_obj, group.by = target_res_para)+
  scale_color_manual(
    values = data_bar_plot_new
  )+
  theme(
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank()
  )
p2 = plot_cells(cds, color_cells_by = "pseudotime",
                label_groups_by_cluster = F,
                split_cells_by = "patient",
                label_branch_points = F,
                label_leaves = F,
                label_roots = F,
                show_trajectory_graph = F,
                group_label_size = 5
)+
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank()
  )+
  scale_colour_gradient2(low = "grey88", mid = "grey88", high = "magenta", midpoint = 15)+
  guides(color=guide_colorbar(title="pseudotime"))

# Add pseudotime to metadata
ccRCC_obj@meta.data$pseudotime <- cds@principal_graph_aux$UMAP$pseudotime[rownames(ccRCC_obj@meta.data)]

# Discretize the pseduotime into 3 pseduotime states (Early, Intermediate and Late)
ccRCC_obj@meta.data$pseudotime_levels <- 
  ifelse(ccRCC_obj@meta.data$pseudotime<quantile(ccRCC_obj@meta.data$pseudotime, 0.33), "Early time", "Intermediate")
ccRCC_obj@meta.data$pseudotime_levels <- 
  ifelse(ccRCC_obj@meta.data$pseudotime>quantile(ccRCC_obj@meta.data$pseudotime, 0.66), "Late time", ccRCC_obj@meta.data$pseudotime_levels)

# Calculating Chi square statitics between discretized pseudotime and tumor stage
ChiRes = chisq.test(ccRCC_obj@meta.data$stage, ccRCC_obj@meta.data$pseudotime_levels)
corrplot_mat = ChiRes$residuals
corrplot_df = as.data.frame(corrplot_mat)
names(corrplot_df) = c("Stages", "Pseudotime Levels","Residuals")
corrplot_df$Stages <- factor(corrplot_df$Stages, levels = c("S1", "S3", "Me"))
corrplot_df$`Pseudotime Levels` <- factor(corrplot_df$`Pseudotime Levels`, 
                                          levels = c("Early time", "Intermediate", "Late time"))
anno_text = ifelse(ChiRes$p.value<0.05, "P < 0.05", paste0("P = ", round(ChiRes$p.value, 3)))
p3 = ggplot(data = corrplot_df)+
  geom_point(aes(x = Stages, y = `Pseudotime Levels`, size = abs(Residuals), color = Residuals))+
  scale_color_gradient2("Pearson residuals", low = "blue", mid = "grey88", high = "red")+
  scale_size_continuous(range=c(5,15))+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.background = element_blank(), legend.position = "right")+
  annotate("text", x = 1.5, y = 3.5, label = anno_text, size = 5, fontface = "bold")+
  guides(size=guide_legend(""))


###############################################################################
###############################################################################
###############################################################################
# Figure S7C and 7D
###############################################################################
###############################################################################
###############################################################################
###############################################################################

for(current_stage in looping_stages){
  scRNA_obj_stage = subset(scRNA_obj, subset = stage == current_stage)
  target_var = target_res_para
  all_types = 0:max(as.numeric(as.character(scRNA_obj@meta.data[,target_var])))
  all_colors = colorRampPalette(palette_color)(length(all_types))
  all_colors = data_bar_plot_new
  names(all_colors) = all_types
  scRNA_obj_stage@meta.data[,plot_var] = factor(scRNA_obj_stage@meta.data[,target_var], levels = all_types)
  
  # define y title
  if(current_stage==looping_stages[1]){
    subpopulations_y_title = "ccRCC Subclones"
    title_function = eval(parse(text="element_text(size=plot_title_size, face='bold', family = 'sans')"))
  }else{
    title_function = eval(parse(text="element_blank()"))
  }
  subpopulations_umpa = DimPlot(scRNA_obj_stage, reduction = "umap", raster = T, group.by = plot_var)+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[as.character(all_types)])+
    ggtitle(current_stage)+
    ylab(subpopulations_y_title)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
          axis.title.x = element_blank(), axis.title.y = title_function,
          # axis.line = element_line(color = "black", size = axis_line_size,
          #                          arrow = grid::arrow(length = unit(0.5, "cm"), ends = "last")),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=plot_title_size, face="bold", family = "sans"), 
          legend.text=element_text(size=legend_text_size, face="bold", family = "sans"), 
          legend.position= "none", legend.justification = "center")+
    guides(colour = guide_legend(override.aes = list(size=legend_shape_size), 
                                 ncol = 4, byrow = T))
  subpopulations_umpas = append(subpopulations_umpas, list(subpopulations_umpa))
  
  all_types = rev(label_tags)
  all_colors = colorRampPalette(heatmap_umap_colors)(length(all_types))
  names(all_colors) = all_types
  
  # make feature plots
  feature_plot_title = ""
  
  feature_plot = lapply(target_signature_list_names, FUN = function(signature){
    if(current_stage==looping_stages[1]){
      feature_plot_y_title = str_replace_all(signature, "_", " ")
      title_function = eval(parse(text="element_text(size=plot_title_size, face='bold', family = 'sans')"))
    }else{
      feature_plot_y_title = ""
      title_function = eval(parse(text="element_blank()"))
    }
    plot_title = ""
    if(signature == target_signature_list_names[1]){
      plot_title = current_stage
    }
    p = DimPlot(scRNA_obj_stage, reduction = "umap", raster = T,
                group.by = paste0(signature, module_score_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      ylab(feature_plot_y_title)+
      ggtitle(plot_title)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
            axis.title.x = element_blank(),
            axis.title.y = title_function,
            # axis.line = element_line(color = "black", size = axis_line_size,
            #                          arrow = grid::arrow(length = unit(0.5, "cm"), ends = "last")),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size=plot_title_size, face="bold", family = "sans"), 
            #plot.title = element_blank(),
            legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
            legend.justification = "center"
      )+
      guides(colour = guide_legend(override.aes = list(size=legend_shape_size)))
    return(p)
  })
  feature_plots = append(feature_plots, feature_plot)
}

mat = 1:(length(looping_stages)*length(target_signature_list_names))
mat = matrix(mat, nrow = length(target_signature_list_names), ncol = length(looping_stages))
feature_plots = feature_plots[c(t(mat))]

subpopulations_umpas_panel = ggarrange(
  plotlist = subpopulations_umpas,
  nrow = 1,
  ncol = length(looping_stages),
  common.legend = T,
  legend = "none"
)
subpopulations_legned = get_legend(subpopulations_umpas)
subpopulations_legned_plot = as_ggplot(subpopulations_legned)
# subpopulations_umpas[[(length(looping_stages)+1)]] = subpopulations_legned_plot
feature_umpas_panel = ggarrange(
  plotlist = feature_plots,
  nrow = length(target_signature_list_names),
  ncol = length(looping_stages),
  common.legend = T,
  legend = "bottom"
)

###############################################################################
###############################################################################
###############################################################################
# Figure 7D Violin quantifications
###############################################################################
###############################################################################
###############################################################################
###############################################################################
pvalue_text_size = 4
stage_colors_dic = c("yellow","orange", "red", "magenta")
names(stage_colors_dic) = stage_dic
stage_colors_template = as.vector(stage_colors_dic[looping_stages])
library(ggpubr)
library(ggprism)
all_vln_plots = lapply(signature_list_umap_names, function(x){
  if(x == signature_list_umap_names[length(signature_list_umap_names)]){
    axis_x_text_function = eval(parse(text="element_text(size=axis_text_size, face = 'bold', family = 'sans')"))
  }else{
    axis_x_text_function = eval(parse(text="element_blank()"))
  }
  # plot_title =str_replace_all(x, "_", " ")
  plot_title = ""
  upper_base = max(scRNA_obj@meta.data[,paste0(x, module_score_tag)])*(0.7)
  lower_base = min(scRNA_obj@meta.data[,paste0(x, module_score_tag)])
  interval = 0.15
  p_table = compare_means(as.formula(paste0(paste0(x, module_score_tag), "~", "stage")), 
                          data = scRNA_obj@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table = p_table[-2,]
  p_table$y.position = seq(upper_base,upper_base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = tibble::as_tibble(p_table)
  reordered_stage = factor(scRNA_obj$stage, levels = looping_stages)
  score = scRNA_obj@meta.data[,paste0(x, module_score_tag)]
  p = ggplot(scRNA_obj@meta.data, aes(x=reordered_stage, y=score))+
    geom_boxplot(aes(fill=reordered_stage), outlier.shape = NA)+
    scale_fill_manual(values = stage_colors_template)+
    # add_pvalue(p_table, tip.length = 0, label.size=pvalue_text_size)+
    ylim(lower_base,upper_base+(interval*2))+
    scale_y_continuous(breaks=c(0,0.5,1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
          axis.line = element_line(color = "black", size = axis_line_size,
                                   arrow = grid::arrow(length = unit(0.1, "cm"), ends = "last")),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = axis_x_text_function, 
          axis.text.y = element_text(size=axis_text_size, face = "bold", family = "sans"), 
          legend.position = "none",
          legend.title = element_blank(),
          legend.key.height = unit(0.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'),
          legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
          legend.box.margin = unit(c(0,0,0,0), "cm"),
          legend.margin = margin(t = 1, b = 1),
          # plot.title = element_text(color = "black", size = plot_title_size, face = "bold", hjust = 0.5, family = "sans")
          plot.title = element_blank()
    )+
    ylab(paste0(x, module_score_tag))
  # guides(fill = guide_legend(override.aes = list(size=legend_shape_size)))+
  p = patchwork::wrap_plots(p)
  return(p)
})
violins = ggarrange(
  plotlist = all_vln_plots,
  ncol = 1,
  nrow = length(signature_list_umap_names),
  common.legend = T,
  legend = "none"
)

###############################################################################
###############################################################################
###############################################################################
# Figure S7D
###############################################################################
###############################################################################
###############################################################################
###############################################################################
scRNA_obj_resampled_meta_df_list = readRDS(
  paste0(intermediate_dir, "columbia_british_individually_harmony_ccRCC_resampled_harmony_meta_df_100sample.rds") 
)

target_res=0.53
target_res_para = paste0("RNA_snn_res.", target_res)
JaccardMatrixMaxTotal = c()
for(i in seq(1,100)){
  scRNA_obj_resampled_meta_df = scRNA_obj_resampled_meta_df_list[[i]]
  resampled_cluster = scRNA_obj_resampled_meta_df[,target_res_para]
  names(resampled_cluster) = rownames(scRNA_obj_resampled_meta_df)
  original_cluster = scRNA_obj@meta.data[,target_res_para]
  names(original_cluster) = rownames(scRNA_obj@meta.data)
  JaccardMatrix = PairWiseJaccardSets(resampled_cluster, original_cluster)
  JaccardMatrixMax = unlist(lapply(data.frame(JaccardMatrix), max))
  JaccardMatrixMaxTotal = rbind(JaccardMatrixMaxTotal, JaccardMatrixMax)
}

JaccardMatrixMaxTotalDf = data.frame(JaccardMatrixMaxTotal, check.names = F)
each_res_mean_of_median = mean(unlist(lapply(JaccardMatrixMaxTotalDf, median)))
rownames(JaccardMatrixMaxTotalDf) = 1:nrow(JaccardMatrixMaxTotalDf)
colnames(JaccardMatrixMaxTotalDf) = seq(0,(ncol(JaccardMatrixMaxTotalDf)-1))
JaccardMatrixMaxTotalDf$id = rownames(JaccardMatrixMaxTotalDf)
JaccardMatrixMaxTotalDfWide = tidyr::gather(JaccardMatrixMaxTotalDf, "vars", "values", -id)

all_types = paste0("sc", 0:14)
all_colors = data_bar_plot_new
names(all_colors) = all_types
JaccardMatrixMaxTotalDfWide$vars = paste0("sc", JaccardMatrixMaxTotalDfWide$vars)
JaccardMatrixMaxTotalDfWide$vars = factor(JaccardMatrixMaxTotalDfWide$vars, levels = all_types)
JaccardBoxPlot = ggplot(JaccardMatrixMaxTotalDfWide)+
  geom_boxplot(mapping = aes(x=vars,y=values,fill = vars))+
  scale_fill_manual(
    breaks = all_types, 
    values=all_colors[as.character(all_types)])+
  xlab("Clusters")+
  ylab("Max Jaccard Index")+
  theme_light()+
  theme(
    axis.text.x = element_text(size = axis_text_size, face = "bold", family = "sans"),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size = axis_text_size, face = "bold", family = "sans"),
    # axis.title.x = element_text(size = axis_title_size, face = "bold", family = "sans"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = axis_title_size, face = "bold", family = "sans"),
    plot.title = element_text(size=plot_title_size, face="bold", hjust = 0.5, family = "sans"),
    legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
    legend.title = element_blank(),
    legend.position = "none",
    legend.box.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(t = 1, b = 1, r = 0, l = 0)
  )+
  scale_x_discrete(labels = c("sc: 0", as.character(1:14)))+
  guides(fill = guide_legend(override.aes = list(size=legend_shape_size), nrow = 2, byrow = T))

###############################################################################
###############################################################################
###############################################################################
# Figure 7E
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# extract pcc df
pcc_expression_df = scRNA_obj@meta.data[,paste0(signature_list_names, module_score_tag)]
bars = lapply(pcc_expression_df, mean)

pcc_revised = function(xv, yv, oxv, oyv, xbar=NULL, ybar=NULL){
  if(is.null(xbar)){
    xbar = mean(xv)
    ybar = mean(yv)
  }
  cov = sum((xv-xbar)*(yv-ybar))
  std_product = sum((oxv-xbar)^2)*sum((oyv-ybar)^2)
  return(cov/sqrt(std_product))
}
combns = list(
  c(2,3),
  c(3,4),
  c(4,5),
  c(5,3)
)
total_res = c()
for(combn in combns){
  each_column = c()
  for(cluster_idx in sort(as.numeric(as.character(unique(scRNA_obj@meta.data[,target_res_para]))))){
    cells.use = rownames(scRNA_obj@meta.data[scRNA_obj@meta.data[,target_res_para]==cluster_idx, ])
    pcc_expression_df_sub = scRNA_obj@meta.data[cells.use, paste0(signature_list_names, module_score_tag)]
    pcc_r = pcc_revised(pcc_expression_df_sub[,names(bars)[combn[1]]], pcc_expression_df_sub[,names(bars)[combn[2]]],
                        pcc_expression_df[,names(bars)[combn[1]]], pcc_expression_df[,names(bars)[combn[2]]], 
                        bars[[combn[[1]]]], bars[[combn[[2]]]])
    each_column[length(each_column)+1] = pcc_r
  }
  total_res = cbind(total_res, each_column)
}

total_res = data.frame(total_res)
df_names = sapply(combns, function(x){
  s1 = signature_list_names[x[1]]
  s2 = signature_list_names[x[2]]
  paste0(s1, "_vs._", s2)
})
names(total_res) = df_names
rownames(total_res) = 0:14
total_res_scaled = total_res

# draw pcc residuals onto UMAPs
# left join meta data with pcc results
total_res_scaled$clusters = as.factor(rownames(total_res_scaled))
join_key= "clusters"
names(join_key) = target_res_para
# scRNA_obj@meta.data = left_join(scRNA_obj@meta.data, total_res_scaled, by=setNames(nm=target_res_para, "clusters"))
# reset PCC variables if they are already there
if(any(names(total_res_scaled)[1:length(df_names)] %in% names(scRNA_obj@meta.data))){
  scRNA_obj@meta.data = scRNA_obj@meta.data[,!(names(scRNA_obj@meta.data) %in% df_names)]
}
meta = dplyr::left_join(scRNA_obj@meta.data, total_res_scaled, by=join_key)
total_res_scaled = total_res_scaled[,-ncol(total_res_scaled)]
scRNA_obj@meta.data[,names(total_res_scaled)] = meta[,names(total_res_scaled)]


all_ccs = lapply(total_res, sum)


pcc_umaps = lapply(names(total_res_scaled), FUN = function(x){
  cor_value = as.numeric(scRNA_obj@meta.data[,x])
  high_cut = sort(unique(cor_value))[length(unique(cor_value))-3]
  if(min(cor_value)<0){
    rescale_values = scales::rescale(c(min(cor_value), 0, high_cut, max(cor_value)))
    rescale_colors = c("blue", "gray88","red", "red")
  }else{
    rescale_values = scales::rescale(c(min(cor_value), high_cut, max(cor_value)))
    rescale_colors = c("gray88","red", "red")
  }
  p = FeaturePlot(scRNA_obj, features = x, label = F, raster = T)+
    # scale_color_gradientn(colours = rescale_colors, values = rescale_values)+
    scale_colour_gradient2(low = "blue", mid = "gray88", high = "red", midpoint = 0,
                           limits = c(-0.04, 0.04), oob = scales::squish)+
    ggtitle(paste0(str_replace_all(x, "_", " "), "\n", "CC = ", round(all_ccs[[x]], 2)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          # axis.line = element_line(color = "black", size = axis_line_size,
          #                          arrow = grid::arrow(length = unit(0.5, "cm"), ends = "last")),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=plot_title_size, face="bold", family = "sans"),
          legend.text = element_text(size=legend_text_size, face="bold", family = "sans"),
          legend.position= "bottom", legend.justification = "center",
          legend.key.width = unit(1.1,"cm"),
          legend.key.height = unit(0.35,"cm"),
    )
  # p = Seurat::LabelClusters(p, id = "ident", box = T, size = umap_label_size, repel = T, force = 10)
  # p = p+scale_fill_manual(values = rep("white",length(unique(Idents(scRNA_obj)))))
  return(p)
})

# output all plots here
plot_dir = "./plots/"
check_dir(plot_dir)
each_width_ratio = 1
png(paste0(plot_dir, "Figure7_A_B_C.png"), width = 15, height = 5, units = "in", res = 300)
ggarrange(
  p1,p2,p3,
  ncol = 3,
  widths = c(2,2.2,2.5)
)
dev.off()

png(paste0(plot_dir, "Figure7_D.png"), width = 8, height = 8, units = "in", res = 300)
ggarrange(
  feature_umpas_panel,
  violins,
  nrow = 1,
  ncol = 2,
  widths = c(length(looping_stages)*each_width_ratio, 0.8)
)
dev.off()

pcc_umaps_panel = ggarrange(
  plotlist = pcc_umaps,
  ncol = 1,
  nrow = length(pcc_umaps),
  common.legend = T,
  legend = "bottom"
)
png(paste0(plot_dir, "Figure7_E.png"),width = 2.6, height = 8, units = "in", res = 300)
pcc_umaps_panel
dev.off()


###############################################################################
###############################################################################
###############################################################################
# Code for Supplementary Figures
###############################################################################
###############################################################################
###############################################################################
###############################################################################

batch_correction_method = batch_correction_methods[3]

# Define cohort names from 1 to 2
cohort_names = paste0("Cohort ", seq(1,2,1))
# define the best resolution for each cluster
res_set = c(0.05, 0.05)
cancer_cell_clusters = c(0,4)

DR_clustering_plots = list()
patient_plots = list()
CA9_plots = list()

for(i in 1:2){
  study_source = study_sources[i]
  # read processed tumor sample data obj
  scRNA_obj = readRDS(paste0("../scRNA_data/",paste(study_source, batch_correction_method, n_of_pc, 
                                                    "PCs", "azimuthed.rds", sep="_")))
  if(i == 1){
    DR_clustering_plot_title = "Clusters"
    patient_plot_title = "Patient"
    CA9_plot_title = "CA9"
  }else{
    DR_clustering_plot_title = ""
    patient_plot_title = ""
    CA9_plot_title = ""
  }
  
  ccRCC_cells_cluster = c(0,4)
  target_res = res_set[i]
  target_res_para = paste0("RNA_snn_res.", target_res)
  group_by_para = "target"
  #scRNA_obj@meta.data[,group_by_para] = factor(scRNA_obj@meta.data[,target_res_para], levels = all_types)
  scRNA_obj@meta.data[,group_by_para] = factor(ifelse(scRNA_obj@meta.data[,target_res_para]==ccRCC_cells_cluster[i],
                                                      "ccRCC cells",
                                                      "Other cells"), levels = c("ccRCC cells", "Other cells"))
  
  all_types = sort(unique(scRNA_obj@meta.data[,group_by_para]))
  all_colors = c("red", "gray")
  names(all_colors) = all_types
  
  DR_clustering_plot = DimPlot(scRNA_obj, reduction = "umap", raster = T, group.by = group_by_para)+
    scale_color_manual(
      breaks = all_types,
      values=all_colors[as.character(all_types)])+
    ylab(cohort_names[i])+
    ggtitle(DR_clustering_plot_title)+
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
  DR_clustering_plots[[length(DR_clustering_plots)+1]] = DR_clustering_plot
  
  # dim plot patient
  plot_var = "patient"
  all_types = unique(scRNA_obj@meta.data[,plot_var])
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  patient_plot = DimPlot(scRNA_obj, reduction = "umap", group.by = plot_var, raster = T)+
    scale_color_manual(
      breaks = all_types,
      values=all_colors[as.character(all_types)])+
    ggtitle(patient_plot_title)
  patient_plots[[length(patient_plots)+1]] = patient_plot
  
  
  target_gene = "CA9"
  DefaultAssay(scRNA_obj) = "RNA"
  scRNA_obj = NormalizeData(scRNA_obj)
  CA9_plot = FeaturePlot(scRNA_obj, features = target_gene, raster = T)+
    scale_color_gradientn(colors = c("grey", "red", "red"), values = c(0,0.4,1))+
    ggtitle(CA9_plot_title)
  CA9_plots[[length(CA9_plots)+1]] = CA9_plot
}

# Run InferCNV
plan("multicore", workers = number_of_thread)
reference_pt = readRDS(
  paste0(intermediate_dir,"reference_pt_cells.rds")
)
#define grouping variable
scRNA_obj$individual_anno = as.character(scRNA_obj$RNA_snn_res.0.53)
grouping_var = "individual_anno"
reference_pt@meta.data[,grouping_var] = "Proximal Tubule"

# restricted to columbia cohort
cells.use <- rownames(scRNA_obj@meta.data[scRNA_obj@meta.data$study=="columbia",])
scRNA_obj <- subset(scRNA_obj, cells = cells.use)

# Define your inferCNV input and output dir
infercnv_input = ""
infercnv_output = ""
# check if infercnv dir exists 
if (!dir.exists(infercnv_input)){
  dir.create(infercnv_input, recursive = T) 
}
if (!dir.exists(infercnv_output)){
  dir.create(infercnv_output, recursive = T) 
}

obj = merge(scRNA_obj, reference_pt)
InfercnvInputs = prepare_infercnv_data(obj, annotation = grouping_var, infercnv_dir = infercnv_input)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=InfercnvInputs[[1]],
                                    annotations_file=InfercnvInputs[[2]],
                                    delim="\t",
                                    gene_order_file=InfercnvInputs[[3]],
                                    ref_group_names = c("Proximal Tubule")
)


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=infercnv_output, 
                             useRaster = F,
                             #k_obs_groups = 14,
                             cluster_by_groups=T, denoise=TRUE,
                             HMM = T, num_threads = number_of_thread)




png(paste0(plot_dir, "FigureS7_A.png"), width = 11, height = 5,
    res = 300, units = "in")
ggarrange(
  DR_clustering_plots,
  CA9_plots,
  patient_plots,
  ncol = 3,
  nrow = 1,
  widths = c(1,1.25,1)
)
dev.off()

# FigureS7_B is the "intercnv.png" from infercnv running results 

png(paste0(plot_dir, "FigureS7_C.png"), width = 6, height = 2, units = "in", res = 300)
ggarrange(plotlist = subpopulations_umpas, legend = "none",
          nrow = 1, ncol = (length(looping_stages)), 
          widths = rep(each_width_ratio, length(looping_stages)))
dev.off()

png(paste0(plot_dir, "FigureS7_D.png"),width = 7, height = 3, units = "in", res = 300)
JaccardBoxPlot
dev.off()

# Individual complexes umaps 
signature_list = list(
  "Complex_I" = signature_list_master$Complex_I,
  "Complex_II" = signature_list_master$Complex_II,
  "Complex_III" = signature_list_master$Complex_III,
  "Complex_IV" = signature_list_master$Complex_IV,
  "Complex_V" = signature_list_master$Complex_V,
  "HIF" = signature_list_master$HIF1A,
  "NRF2" = signature_list_master$NRF2
)

# convert gene name from lower case to upper
signature_list = lapply(signature_list, toupper)
# extract names for this signature list
signature_list_names = names(signature_list)

# to profile the expression pattern of pathway signatures,
# we use the original RNA assay
DefaultAssay(scRNA_obj) = "RNA"
scRNA_obj = NormalizeData(scRNA_obj)


# Add module score for for all pathway signatures
scRNA_obj = AddModuleScore(scRNA_obj,
                           features = signature_list,
                           assay = assay_in_use,
                           nbin = 24,
                           ctrl = 100,
                           name = signature_list_names)
scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)] =
  scRNA_obj@meta.data[, paste0(signature_list_names, 1:length(signature_list_names))]

# Add z-score transformed module score
if(length(signature_list_names)==1){
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")] = 
    scale(scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)])
}else{
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")] =
    as.data.frame(lapply(scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag)], scale))
}

#remove redundant features
remove_column_index = match(paste0(signature_list_names, 1:length(signature_list_names)),
                            names(scRNA_obj@meta.data))
scRNA_obj@meta.data =
  scRNA_obj@meta.data[, -remove_column_index]


break_cutoffs = c(-Inf,-1.5,-0.5, 0.5, 1.5, Inf)
label_tags = create_label_tag(break_cutoffs)
if (length(signature_list_names) == 1) {
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_group")] =
    sapply(
      scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")],
      FUN = function(vx) {
        vx = cut(vx, breaks = break_cutoffs, labels = label_tags)
      }
    )
} else{
  scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_group")] =
    lapply(
      scRNA_obj@meta.data[, paste0(signature_list_names, module_score_tag, "_zscore")],
      FUN = function(vx) {
        vx = cut(vx, breaks = break_cutoffs, labels = label_tags)
      }
    )
}


# loop each stage and get umaps individually
plot_title_size = 15
axis_title_size = 15
legend_text_size = 12
legend_shape_size = 5
axis_line_size = 0.5
axis_text_size = 12
umap_label_size = 1.5
target_signature_list_names = signature_list_names
subpopulations_umpas = list()
feature_plots = list()
#define the dummy plot var to avoid changing of original data
plot_var = "target"
pmt = 0
pmr = 0
pmb = 0
pml = 0
for(current_stage in looping_stages){
  scRNA_obj_stage = subset(scRNA_obj, subset = stage == current_stage)
  target_var = target_res_para
  all_types = 0:max(as.numeric(as.character(scRNA_obj@meta.data[,target_var])))
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  scRNA_obj_stage@meta.data[,plot_var] = factor(scRNA_obj_stage@meta.data[,target_var], levels = all_types)
  
  # define y title
  if(current_stage==looping_stages[1]){
    subpopulations_y_title = "ccRCC Subclones"
    title_function = eval(parse(text="element_text(size=plot_title_size, face='bold', family = 'sans')"))
  }else{
    title_function = eval(parse(text="element_blank()"))
  }
  subpopulations_umpa = DimPlot(scRNA_obj_stage, reduction = "umap", raster = T, group.by = plot_var)+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[as.character(all_types)])+
    ggtitle(current_stage)+
    ylab(subpopulations_y_title)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
          axis.title.x = element_blank(), axis.title.y = title_function,
          # axis.line = element_line(color = "black", size = axis_line_size,
          #                          arrow = grid::arrow(length = unit(0.5, "cm"), ends = "last")),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=plot_title_size, face="bold", family = "sans"), 
          legend.text=element_text(size=legend_text_size, face="bold", family = "sans"), 
          legend.position= "none", legend.justification = "center")+
    guides(colour = guide_legend(override.aes = list(size=legend_shape_size), 
                                 ncol = 4, byrow = T))
  
  subpopulations_umpas = append(subpopulations_umpas, list(subpopulations_umpa))
  
  all_types = rev(label_tags)
  all_colors = colorRampPalette(heatmap_umap_colors)(length(all_types))
  names(all_colors) = all_types
  
  # make feature plots
  
  feature_plot = lapply(target_signature_list_names, FUN = function(signature){
    if(current_stage==looping_stages[1]){
      feature_plot_y_title = str_replace_all(signature, "_", " ")
      title_function = eval(parse(text="element_text(size=plot_title_size, face='bold', family = 'sans')"))
    }else{
      feature_plot_y_title = ""
      title_function = eval(parse(text="element_blank()"))
    }
    plot_title = ""
    if(signature == target_signature_list_names[1]){
      plot_title = current_stage
    }
    p = DimPlot(scRNA_obj_stage, reduction = "umap", raster = T,
                group.by = paste0(signature, module_score_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      ylab(feature_plot_y_title)+
      ggtitle(plot_title)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
            axis.title.x = element_blank(),
            axis.title.y = title_function,
            # axis.line = element_line(color = "black", size = axis_line_size,
            #                          arrow = grid::arrow(length = unit(0.5, "cm"), ends = "last")),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size=plot_title_size, face="bold", family = "sans"), 
            #plot.title = element_blank(),
            legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
            legend.justification = "center"
      )+
      guides(colour = guide_legend(override.aes = list(size=legend_shape_size)))
    return(p)
  })
  feature_plots = append(feature_plots, feature_plot)
}

mat = 1:(length(looping_stages)*length(target_signature_list_names))
mat = matrix(mat, nrow = length(target_signature_list_names), ncol = length(looping_stages))
feature_plots = feature_plots[c(t(mat))]

subpopulations_umpas_panel = ggarrange(
  plotlist = subpopulations_umpas,
  nrow = 1,
  ncol = length(looping_stages),
  common.legend = T,
  legend = "none"
)
subpopulations_legned = get_legend(subpopulations_umpas)
subpopulations_legned_plot = as_ggplot(subpopulations_legned)
# subpopulations_umpas[[(length(looping_stages)+1)]] = subpopulations_legned_plot
feature_umpas_panel = ggarrange(
  plotlist = feature_plots,
  nrow = length(target_signature_list_names),
  ncol = length(looping_stages),
  common.legend = T,
  legend = "none"
)


# plot out violin

pvalue_text_size = 4
#legend_text_size = 7
stage_colors_dic = c("yellow","orange", "red", "magenta")
names(stage_colors_dic) = stage_dic
stage_colors_template = as.vector(stage_colors_dic[looping_stages])
library(ggpubr)
library(ggprism)
all_vln_plots = lapply(signature_list_names, function(x){
  if(x == signature_list_names[length(signature_list_names)]){
    axis_x_text_function = eval(parse(text="element_text(size=axis_text_size, face = 'bold', family = 'sans')"))
  }else{
    axis_x_text_function = eval(parse(text="element_blank()"))
  }
  # plot_title =str_replace_all(x, "_", " ")
  plot_title = ""
  upper_base = max(scRNA_obj@meta.data[,paste0(x, module_score_tag)])*(0.7)
  scale_break = c(0,0.5,1,1.5)
  interval = 0.2
  if(x == "Complex_II"){
    upper_base = upper_base-0.5
    scale_break = c(0,0.5,1)
    interval = 0.3
  }
  if(x == "Complex_III"){
    interval = 0.3
  }
  if(x=="Complex_V"){
    interval = 0.3
  }
  lower_base = min(scRNA_obj@meta.data[,paste0(x, module_score_tag)])
  p_table = compare_means(as.formula(paste0(paste0(x, module_score_tag), "~", "stage")), 
                          data = scRNA_obj@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table = p_table[-2,]
  p_table$y.position = seq(upper_base,upper_base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = tibble::as_tibble(p_table)
  reordered_stage = factor(scRNA_obj$stage, levels = looping_stages)
  score = scRNA_obj@meta.data[,paste0(x, module_score_tag)]
  p = ggplot(scRNA_obj@meta.data, aes(x=reordered_stage, y=score))+
    geom_boxplot(aes(fill=reordered_stage), outlier.shape = NA)+
    scale_fill_manual(values = stage_colors_template)+
    add_pvalue(p_table, tip.length = 0, label.size=pvalue_text_size)+
    ylim(lower_base,upper_base+(interval*2))+
    scale_y_continuous(breaks=scale_break)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.margin = unit(c(pmt,pmr,pmb, pml), "cm"),
          axis.line = element_line(color = "black", size = axis_line_size,
                                   arrow = grid::arrow(length = unit(0.1, "cm"), ends = "last")),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = axis_x_text_function, 
          axis.text.y = element_text(size=axis_text_size, face = "bold", family = "sans"), 
          legend.position = "none",
          legend.title = element_blank(),
          legend.key.height = unit(0.4, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'),
          legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
          legend.box.margin = unit(c(0,0,0,0), "cm"),
          legend.margin = margin(t = 1, b = 1),
          # plot.title = element_text(color = "black", size = plot_title_size, face = "bold", hjust = 0.5, family = "sans")
          plot.title = element_blank()
    )+
    ylab(paste0(x, module_score_tag))
  # guides(fill = guide_legend(override.aes = list(size=legend_shape_size)))+
  p = patchwork::wrap_plots(p)
  return(p)
})

violins = ggarrange(
  plotlist = all_vln_plots,
  ncol = 1,
  nrow = length(signature_list_names),
  common.legend = T,
  legend = "none"
)

Individual_complexes_panel = ggarrange(
  feature_umpas_panel,
  violins,
  nrow = 1,
  ncol = 2,
  widths = c(length(looping_stages), 0.8)
)

png(paste0(plot_dir, "FigureS7_E.png"), width = 8, height = 12, units = "in", res = 300)
Individual_complexes_panel
dev.off()






