obj <- readRDS("./umap_tunned_103530_cells.rds")
library(Seurat)
library(ggplot2)
set.seed(42)

# UMAP降维查看批次效应
obj <- RunUMAP(
  obj,
  dims = 1:12,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.2,
  #min_dist must be less than or equal to spread
  reduction = "integrated.jointpca",
  reduction.name = "jointpca.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(obj, group.by = "Batch", reduction = "jointpca.umap") + ggtitle("integrated by jointpca")

obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 200,
  spread = 0.25,
  min.dist = 0.25,
  #min_dist must be less than or equal to spread
  reduction = "integrated.harmony",
  reduction.name = "harmony.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(obj, group.by = "Batch", reduction = "harmony.umap") + ggtitle("integrated by harmony")

# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:12,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.25,
  #min_dist must be less than or equal to spread
  reduction = "integrated.rpca",
  reduction.name = "rpca.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(obj, group.by = "Batch", reduction = "rpca.umap") + ggtitle("integrated by RPCA")

# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 100,
  spread = 0.35,
  min.dist = 0.35,
  #min_dist must be less than or equal to spread
  reduction = "integrated.cca",
  reduction.name = "cca.umap",
  umap.method = "umap-learn",
  metric = "cosine"
)

DimPlot(obj, group.by = "Batch", reduction = "cca.umap") + ggtitle("integrated by CCA")

# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:30,
  # n.neighbors = 100,
  # spread = 0.25,
  # min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "integrated.scvi",
  reduction.name = "scvi.umap",
  # umap.method = "umap-learn",
  metric = "cosine"
)

DimPlot(obj, group.by = "Batch", reduction = "scvi.umap") + ggtitle("integrated by SCVI")

# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:12,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.25,
  #min_dist must be less than or equal to spread
  reduction = "pca",
  reduction.name = "pca.umap",
  umap.method = "umap-learn",
  metric = "cosine"
)

DimPlot(obj, group.by = "Batch", reduction = "pca.umap") + ggtitle("unintegrated")

# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:12,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.25,
  #min_dist must be less than or equal to spread
  reduction = "integrated.mnn",
  reduction.name = "mnn.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(obj, group.by = "Batch", reduction = "mnn.umap") + ggtitle("mnn")

saveRDS(obj, "./umap_tunned_103530_cells.rds")


library(Seurat)
library(ggplot2)
library(cowplot)  # 或者使用 patchwork: library(patchwork)

# 生成每个 UMAP 的图形
p1 <- DimPlot(obj,
              reduction = "mnn.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("MNN UMAP")
p2 <- DimPlot(obj,
              reduction = "jointpca.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("JointPCA UMAP")
p3 <- DimPlot(obj,
              reduction = "pca.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("PCA UMAP")
p4 <- DimPlot(obj,
              reduction = "harmony.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("Harmony UMAP")
p5 <- DimPlot(obj,
              reduction = "rpca.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("RPCA UMAP")
p6 <- DimPlot(obj,
              reduction = "cca.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("CCA UMAP")
p7 <- DimPlot(obj,
              reduction = "scvi.umap",
              group.by = "jointpca_cell_type",
              raster = F) + ggtitle("scVI UMAP")

# 将多个图组合在一起
combined_plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, ncol = 3)

# 保存图像到文件
ggsave(
  "combined_umap_plots.png",
  plot = combined_plot,
  width = 16,
  height = 12
)
