# ccRCC scRNA-seq analysis
## Enviroment setups
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
- scclusteval
## Usage
* Run scripts using processed data (**Recommended**): If you do not want to start from the very raw data,please use [this link](https://drive.google.com/file/d/18sO6wuIoghZ1LsVL46PgMi_Sj4Amy67f/view?usp=sharing) and [this link](https://drive.google.com/file/d/1eTA7yCiLVzVfXBMDxtyEDRtNsK_uANYx/view?usp=sharing) to download the data to "intermediate_files" dir and run "Manuscript_figures.R".

* Run scripts using raw data: If you want to start from the very raw data, please use [this link](https://drive.google.com/file/d/1CEdUnspMfthex0XSxemxCy1j5TWkWxCE/view?usp=sharing) to download the "raw_data.rds" and run "Preprocessing.R", "Identify_ccRCC_cells.R", "Intergrate_ccRCC_cells_datasets.R" and "Manuscript_figures.R" accordingly.
## System recommendation
* It is recommended to use Linux-like systems (Linux or MacOs) as the platform and have memory>=30G, if you would like to run scripts using raw data
