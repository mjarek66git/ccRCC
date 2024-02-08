# number of threads to use
number_of_thread = 4
# mt threshold
mt_threshold = 10
# number of highly variable genes to select
n_of_highly_variable_genes = 3000
# batch correction methods
batch_correction_methods = c("SeuratSCT", "Seurat", "harmony", "scVI")
# study sources
study_sources = c("columbia", "british", "harvard", "michigan")
# number of PCA space to select for downstream analysis
n_of_pc = 30
# clustering resolution scope
start_res = 0.01
end_res = 1
step_res = 0.01
resolution_scope = seq(start_res, end_res, step_res)

#color settings
palette_color = c("firebrick1", "darkorange", "gold", 
                  "chartreuse", "blue", "mediumorchid1", 
                  "cyan", "lightpink", "sienna4")
heatmap_umap_colors = c("gold", "gold3", "grey", "cyan", "blue")
# Define the directory path you want to create
intermediate_dir <- "./intermediate_files/"


