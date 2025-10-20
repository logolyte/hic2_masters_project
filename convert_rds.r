# Analyze the rds file

library(Seurat)
library(SeuratDisk)

object <- readRDS("/home/love/Projects/hic2_masters_project/RNA Sequencing Data/220701Pool1-Pool3.rds")
SaveH5Seurat(object, filename = "220701Pool1-Pool3.h5Seurat")
Convert("220701Pool1-Pool3.h5Seurat", dest = "h5ad")

