library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
set.seed(42)

filtered_obj <- readRDS("./filtered_jointpca_annotated_99252_cells.rds")

DimPlot(filtered_obj, reduction = "cca.umap")
DimPlot(filtered_obj, reduction = "mnn.umap")
