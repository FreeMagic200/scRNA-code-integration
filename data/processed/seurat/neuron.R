library(Seurat)
library(ggplot2)
library(scCustomize)
library(mrtree)
library(clustree)
library(dplyr)
library(SeuratWrappers)
set.seed(42)

options(
  future.globals.maxSize = Inf,
  future.seed = TRUE,
  future.seed = T
)
obj <- readRDS("./Neuron.RDS")

DimPlot(obj)

# preprocess
obj <- FindVariableFeatures(obj)

obj <- ScaleData(obj)

obj <- RunPCA(obj)

# Perform JointPCA integration
obj <- IntegrateLayers(
  object = obj,
  method = JointPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.jointpca",
  verbose = FALSE
)

# Perform scvi integration
obj <- IntegrateLayers(
  object = obj,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/opt/miniforge3/envs/scvi/",
  verbose = FALSE
)

ElbowPlot(obj, ndims = 50)

obj <- RunUMAP(
  obj,
  dims = 1:30,
  # n.neighbors = 15,
  # spread = 0.2,
  # min.dist = 0.2,
  reduction = "integrated.jointpca"
)

DimPlot(obj, reduction = "umap")

FeaturePlot(obj,"Pomc",order = T,raster = F, reduction = "umap")

saveRDS(obj, "./Neuron_integrated.RDS")
gc()

obj <- readRDS("./Neuron_integrated.RDS")

obj <- FindNeighbors(
  obj,
  dims = 1:12,
  reduction = "integrated.jointpca",
  k.param = 5,
  annoy.metric = "euclidean"
)

obj <- FindClusters(
  obj,
  method = "igraph",
  algorithm = 4,
  resolution = seq(0, 2, 0.1),
  verbose = T
)

k5_tree<-clustree(obj,prefix = "RNA_snn_res.")

# 假设你的数据框名为 k5_tree$data
data <- k5_tree$data

# 处理数据
plot_data <- data %>%
  group_by(RNA_snn_res.) %>%
  summarise(
    avg_stability = mean(sc3_stability, na.rm = TRUE),
    cluster_count = n_distinct(cluster)
  ) %>%
  mutate(resolution = as.numeric(as.character(RNA_snn_res.)))

# 创建折线图
ggplot(plot_data, aes(x = resolution, y = avg_stability)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Resolution",
    y = "Average SC3 Stability",
    title = "Average SC3 Stability vs. Resolution"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = plot_data$resolution) +
  geom_text(aes(label = cluster_count), vjust = -0.5, size = 3)

# 保存图片
ggsave("neu_k5_clustree.png", plot = k5_tree, width = 12, height = 14, dpi = 200)

DimPlot(obj, reduction = "umap", group.by = "RNA_snn_res.1.6",raster = F,label = T)

saveRDS(obj,"./Neuron_clustered.RDS")
