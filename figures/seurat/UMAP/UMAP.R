library(Seurat)
library("Nebulosa")
set.seed(42)
obj<-readRDS("../../../data/processed/seurat/06_integrated_data.rds")
obj<- RunUMAP(
  obj,
  dims = 1:30,
  reduction = "integrated.cca",
  n.neighbors = 90,
  min.dist = 0.4,
  spread = 1,
  umap = "umap-learn",
  metric = "correlation",
  reduction.name = "umap.cca"
)

plot_density(obj, "Pdgfra")


