source("tools.R")
source("Hyperparameters_for_pipeline.R")
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(scclusteval)
scRNA_obj = readRDS(
  paste0(intermediate_dir,"columbia_british_individually_harmony_ccRCC_harmony.rds")
)

# p1 = FeaturePlot(scRNA_obj, features = "HIF1A")
# p2 = FeaturePlot(scRNA_obj, features = "NFE2L2")
# p3 = DimPlot(scRNA_obj, group.by = "RNA_snn_res.0.53")
# 
# p1 | p2 |p3
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
# Figure 6A and 6C UMAPs
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
# Figure 6C Violin quantifications
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
# Figure 6B
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
# Figure 6D
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
png(paste0(plot_dir, "Figure6_A.png"), width = 6, height = 2, units = "in", res = 300)
ggarrange(plotlist = subpopulations_umpas, legend = "none",
          nrow = 1, ncol = (length(looping_stages)), 
          widths = rep(each_width_ratio, length(looping_stages)))
dev.off()

png(paste0(plot_dir, "Figure6_B.png"),width = 7, height = 3, units = "in", res = 300)
JaccardBoxPlot
dev.off()

png(paste0(plot_dir, "Figure6_C.png"), width = 8, height = 8, units = "in", res = 300)
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
png(paste0(plot_dir, "Figure6_D.png"),width = 2.6, height = 8, units = "in", res = 300)
pcc_umaps_panel
dev.off()

