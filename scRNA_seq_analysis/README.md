# ccRCC scRNA-seq analysis
## Enviroment setups
### Recommending memory
It is recommended to have memory>=30G to run scripts in this directory
### Required libraries
- Seurat(V4)
- future
- harmony
- Azimuth
- ggplot2
- ggpubr
- dplyr
- stringr
- ggprism
- tibble
- patchwork
- grid
## Usage
If you want to start from the very raw data, please use [this link](https://drive.google.com/file/d/1CEdUnspMfthex0XSxemxCy1j5TWkWxCE/view?usp=sharing) to download the "raw_data.rds" and run "Preprocessing.R", "Identify_ccRCC_cells.R", "Intergrate_ccRCC_cells_datasets.R" and "Manuscript_figures.R" accordingly.
Otherwise, please use [this link](https://drive.google.com/file/d/18sO6wuIoghZ1LsVL46PgMi_Sj4Amy67f/view?usp=sharing) to download the data to "intermediate_files" dir and run "Manuscript_figures.R"