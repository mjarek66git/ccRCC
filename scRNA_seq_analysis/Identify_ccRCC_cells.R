set.seed(10)
library(Seurat)
library(ggplot2)
library(ggpubr)
source("Pipelines.R")
source("Hyperparameters_for_pipeline.R")
source("tools.R")
check_dir(intermediate_dir)
batch_correction_method = batch_correction_methods[3]

# Define cohort names from 1 to 2
cohort_names = paste0("Cohort ", seq(1,2,1))
# define the best resolution for each cluster
res_set = c(0.05, 0.05)
ccRCC_cells_cluster = c(0,4)

DR_clustering_plots = list()
patient_plots = list()
CA9_plots = list()

for(i in 1:2){
  study_source = study_sources[i]
  # read processed tumor sample data obj
  scRNA_obj = readRDS(paste0(intermediate_dir, paste(study_source, batch_correction_method, n_of_pc, 
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
  ## dim plot clustering
  target_res = res_set[i]
  target_res_para = paste0("RNA_snn_res.", target_res)
  group_by_para = "target"
  #scRNA_obj@meta.data[,group_by_para] = factor(scRNA_obj@meta.data[,target_res_para], levels = all_types)
  scRNA_obj@meta.data[,group_by_para] = factor(ifelse(scRNA_obj@meta.data[,target_res_para]==ccRCC_cells_cluster[i],
                                                      "ccRCC cells",
                                                      "Other cells"), levels = c("ccRCC cells", "Other cells"))
  
  saveRDS(scRNA_obj, paste0(intermediate_dir, paste(study_source, batch_correction_method, "ccRCC.rds", sep = "_")))
  
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


legend_figure = get_legend(DR_clustering_plot)
Clustering_legend_ggplot = ggpubr::as_ggplot(legend_figure)

plot_title_size = 15
axis_title_size = 15
umap_label_size = 4
legend_text_size = 12
legend_shape_size = 5
DR_clustering_plots_display = lapply(DR_clustering_plots, function(x){
  # x = Seurat::LabelClusters(x, id = "target", box = T, size = umap_label_size, repel = T, 
  #                           box.padding = 0.15, label.padding = 0.15)+
  #   scale_fill_manual(values = rep("white",length(unique(x$data$target))))
  
  x+theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=axis_title_size, face="bold", family = "sans"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(size=plot_title_size, face="bold", family = "sans"),
    legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
    legend.position= "none",
    legend.box.spacing = unit(0, "cm")
  )
})

patient_plots_display = lapply(patient_plots, function(x){
  x+guides(colour = guide_legend(override.aes = list(size=legend_shape_size), ncol= 1))+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(size=plot_title_size, face="bold", family = "sans"),
      # legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
      legend.text=element_blank(),
      legend.position= "right",
      legend.box.spacing = unit(0, "cm")
    )
})
legend_figures = list()
legend_figure = get_legend(patient_plots_display[[1]])
legend_figures[[length(legend_figures)+1]] = ggpubr::as_ggplot(legend_figure)
legend_figure = get_legend(patient_plots_display[[2]])
legend_figures[[length(legend_figures)+1]] = ggpubr::as_ggplot(legend_figure)

patient_plots_display = lapply(patient_plots_display, function(x){
  x+theme(
    legend.position= "none"
  )
})

CA9_plots_display = lapply(CA9_plots, function(x){
  x+guides(fill = guide_legend(override.aes = list(size=legend_shape_size)))+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(size=plot_title_size, face="bold", family = "sans"),
      legend.text=element_text(size=legend_text_size, face="bold", family = "sans"),
      legend.position= "bottom",
      legend.justification = "center",
      legend.box.spacing = unit(0, "cm")
    )
})

legend_figure = get_legend(CA9_plots_display[[1]])
CA9_legend_ggplot = ggpubr::as_ggplot(legend_figure)

CA9_plots_display = lapply(CA9_plots_display, function(x){
  x+theme(
    legend.position= "none"
  )
})

png("plots/manuscript_plots/Fig7S_Panel_A_1.png", width = 3.3, height = 2.5, units = "in", res = 300)
DR_clustering_plots_display[[1]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_2.png", width = 3.3, height = 2.5, units = "in", res = 300)
DR_clustering_plots_display[[2]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_1_and_2_legend.png", width = 3, height = 0.5, units = "in", res = 300)
Clustering_legend_ggplot
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_3.png", width = 3.3, height = 2.5, units = "in", res = 300)
CA9_plots_display[[1]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_4.png", width = 3.3, height = 2.5, units = "in", res = 300)
CA9_plots_display[[2]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_3_and_4_legend.png", width = 2, height = 0.5, units = "in", res = 300)
CA9_legend_ggplot
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_5.png", width = 3.3, height = 2.5, units = "in", res = 300)
patient_plots_display[[1]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_6.png", width = 3.3, height = 2.5, units = "in", res = 300)
patient_plots_display[[2]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_5_lenged.png", width = 0.5, height = 2, units = "in", res = 300)
legend_figures[[1]]
dev.off()
png("plots/manuscript_plots/Fig7S_Panel_A_6_lenged.png", width = 0.5, height = 2.5, units = "in", res = 300)
legend_figures[[2]]
dev.off()


cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_1.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
DR_clustering_plots_display[[1]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_2.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
DR_clustering_plots_display[[2]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_1_and_2_legend.eps", width = 3, height = 0.5, fallback_resolution = 300)
Clustering_legend_ggplot
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_3.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
CA9_plots_display[[1]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_4.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
CA9_plots_display[[2]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_3_and_4_legend.eps", width = 2, height = 0.5, fallback_resolution = 300)
CA9_legend_ggplot
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_5.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
patient_plots_display[[1]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_6.eps", width = 3.3, height = 2.5, fallback_resolution = 300)
patient_plots_display[[2]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_5_lenged.eps", width = 0.5, height = 2, fallback_resolution = 300)
legend_figures[[1]]
dev.off()
cairo_ps("plots/manuscript_plots/Fig7S_Panel_A_6_lenged.eps", width = 0.5, height = 2.5, fallback_resolution = 300)
legend_figures[[2]]
dev.off()



ccRCC_cells_panel = ggarrange(
  DR_clustering_panel,
  CA9_panel,
  patient_panel,
  ncol = 3,
  nrow = 1,
  widths = c(1,1.25,1)
)
