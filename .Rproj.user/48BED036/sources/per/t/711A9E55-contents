library(Seurat)
library(ggplot2)
library(stringr)
library(SeuratWrappers)
library(future)
library(patchwork)
set.seed(42)
obj <- readRDS("./06_integrated_data.rds")
# obj <- readRDS("./afterscvi.rds")
# DimPlot(obj)
plan("sequential", workers = 1)
options(future.globals.maxSize = Inf, future.seed = TRUE,future.seed = T,)
plan()

# 未整合的 ----------------------------------------------------
# UMAP降维查看批次效应
obj <- RunUMAP(
  obj,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "na.umap"
)
DimPlot(obj, group.by = "Batch", reduction = "na.umap") + ggtitle("not integrated")

# 先查看几个Marker基因
# FeaturePlot(obj, c("Pdgfra", "Mbp"), reduction = "na.umap")

# 聚类
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca",k.param = 30)
obj <- FindClusters(
  obj,
  cluster.name = "pca.clust",
  method = "igraph",
  algorithm = 4,
  resolution = 0.6,
  verbose = T
)

# 可视化聚类结果
DimPlot(obj, reduction = "na.umap", group.by = "pca.clust")

DotPlot(obj, group.by = "pca.clust", features = marker_genes_list) + RotatedAxis() +
  ggtitle("not integrated") +
  theme(plot.title = element_text(hjust = 0.5))

obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 20,
  spread = 0.25,
  min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "pca",
  reduction.name = "na.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(
  obj,
  reduction = "na.umap",
  group.by = "pca.clust",
  raster = F,
  label = T
) + ggtitle("not integrated")


# CCA --------------------------------------------------------
# UMAP降维查看批次效应
obj <- RunUMAP(
  obj,
  dims = 1:22,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "integrated.cca",
  reduction.name = "cca.umap",
  umap.method = "umap-learn",
  metric = "cosine"
)

DimPlot(obj, group.by = "Batch", reduction = "cca.umap") + ggtitle("integrated by CCA")

# 聚类
obj <- FindNeighbors(obj, dims = 1:15, reduction = "integrated.cca",k.param = 60)
obj <- FindClusters(
  obj,
  cluster.name = "cca.clust",
  method = "igraph",
  algorithm = 4,
  resolution = 1.8,
  verbose = T
)

# 可视化聚类结果
DimPlot(obj, reduction = "cca.umap", group.by = "cca.clust")

DotPlot(obj, group.by = "cca.clust", features = marker_genes_list) +
  RotatedAxis() +
  ggtitle("integrated by CCA") +
  theme(plot.title = element_text(hjust = 0.5))

# 创建一个细胞类型到cluster的映射字典 22/18不确定
cca_cell_type_mapping <- list(
  GABA = c(1, 7, 10, 11, 16, 17, 19, 21, 22),
  GLU = c(4, 6, 12, 13, 14, 15, 23, 18),
  IPC_Ascl1 = c(8),
  IPC_Neurog2 = c(20),
  TC_RGC = c(9),
  AS_RGC = c(5),
  RGC = c(2),
  Microglia_Fib = 25,
  OD_OPC = 3,
  EC = 24
)

obj$cca_cell_type <- NA

# 根据映射字典为每个cluster添加cell_types列
for (i in seq_along(cca_cell_type_mapping)) {
  clusters <- cca_cell_type_mapping[[i]]
  cell_type <- names(cca_cell_type_mapping)[i]
  obj$cca_cell_type[obj$cca.clust %in% clusters] <- cell_type
}

DimPlot(
  obj,
  reduction = "cca.umap",
  group.by = "cca_cell_type",
  raster = F,
  label = T
) + ggtitle("integrated by CCA")

obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 30,
  spread = 0.25,
  min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "integrated.cca",
  reduction.name = "cca.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(
  obj,
  reduction = "cca.umap",
  group.by = "cca_cell_type",
  raster = F,
  label = T
) + ggtitle("integrated by CCA")


# RPCA-------------------------------------------------------
# UMAP降维
obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 100,
  spread = 0.25,
  min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "integrated.rpca",
  reduction.name = "rpca.umap",
  umap.method = "umap-learn",
  metric = "cosine"
)

DimPlot(obj, group.by = "Batch", reduction = "rpca.umap") + ggtitle("integrated by RPCA")

# 聚类
# 细分聚类分出OD
obj <- FindNeighbors(obj, dims = 1:15, reduction = "integrated.rpca",k.param = 5,annoy.metric = "cosine")
obj <- FindClusters(
  obj,
  method = "igraph",
  algorithm = 4,
  resolution = seq(2.25,3,0.25),
  verbose = T
)

obj <- FindClusters(
  obj,
  method = "igraph",
  algorithm = 4,
  resolution = c(1.6,1.7,1.8,1.9),
  verbose = T
)

obj <- FindClusters(
  obj,
  method = "igraph",
  algorithm = 4,
  resolution = 1.7,
  verbose = T
)

gc()

# library(clustree)
k5_tree<-clustree(obj,prefix = "RNA_snn_res.")

library(dplyr)
library(ggplot2)

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
ggsave("k5_clustree_indent_1.png", plot = k5_tree, width = 12, height = 14, dpi = 200)

DimPlot(obj, reduction = "rpca.umap", group.by = "RNA_snn_res.1.7",raster = F,label = T)

# # 转换为数值型，并重新设置因子水平
# obj$RNA_snn_res.2.5 <- factor(obj$RNA_snn_res.2.5, levels = as.character(sort(as.numeric(levels(obj$RNA_snn_res.2.5)))))

# 创建点图
p_1.7 <- DotPlot(obj, group.by = "RNA_snn_res.1.7", features = marker_genes_list) + 
  RotatedAxis() +
  ggtitle("integrated by RPCA") +
  theme(plot.title = element_text(hjust = 0.5))

# 保存图形
ggsave("res1.7_RPCA_clust.png", plot = p_1.7, width = 26, height = 12, dpi = 300)

obj$rpca.clust <-obj$RNA_snn_res.1.7

# 转换为数值型，并重新设置因子水平
obj$rpca.clust <- factor(obj$rpca.clust, levels = as.character(sort(as.numeric(levels(obj$rpca.clust)))))

# check quality
plot_density(obj,"nFeature_RNA",reduction = "rpca.umap")
VlnPlot(obj,"nFeature_RNA",group.by = "rpca.clust",pt.size = 0)
plot_density(obj,"percent.mt",reduction = "rpca.umap")

# marker dotplot
rpca_clust_dotplot<-DotPlot(obj, group.by = "rpca.clust", features = unique(unlist(marker_genes_list))) + RotatedAxis() +
  ggtitle("RPCA clust Markers") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

ggsave("rpca_clust_dotplot.png", rpca_clust_dotplot, width = 34, height = 15, dpi = 300)


FeaturePlot(obj,marker_genes_list$AS,order = F,raster = F,max.cutoff = 0.95,pt.size = 0.1)
FeaturePlot(obj,c("Crym","Col25a1","Rax"),order = F,raster = F,max.cutoff = 0.95,pt.size = 0.1)
FeaturePlot(obj,c("Pdgfra"),order = F,raster = F,max.cutoff = 0.95,pt.size = 0.1)

# 获取基因列表
genes_to_plot <- c(marker_genes, "Col25a1")

# 计算基因数量
num_genes <- length(genes_to_plot)

# 动态调整 ncol
ncol <- min(5, ceiling(sqrt(num_genes)))  # 最多5列，但会根据基因数量调整

# 创建 FeaturePlot
p <- FeaturePlot(obj, 
                 features = genes_to_plot, 
                 order = FALSE, 
                 raster = FALSE, 
                 max.cutoff = 0.95, 
                 pt.size = 0.1, 
                 ncol = ncol)
p1 <- FeaturePlot(obj, 
                 features = genes_to_plot, 
                 order = T, 
                 raster = FALSE, 
                 max.cutoff = 0.95, 
                 pt.size = 0.1, 
                 ncol = ncol)

# 计算合理的图片大小
width <- ncol * 5  # 每列5英寸宽
height <- ceiling(num_genes / ncol) * 5  # 每行5英寸高

# 保存图片
ggsave("FeaturePlot_genes.png", plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)
ggsave("FeaturePlot_ordered_genes.png", plot = p1, width = width, height = height, dpi = 300, limitsize = FALSE)

VlnPlot(obj,c("nFeature_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb"),group.by = "rpca.clust",pt.size = 0)

# 创建一个细胞类型到cluster的映射字典 25不确定
rpca_cell_type_mapping <- list(
  GABA = c(3,6,7,10,13,15,16,17,18,22,27),
  GLU = c(1,5,8,9,14,23,25,26,28,29,31,33), #25/26/31很少表达Slc17a6/Slc32a1 25是pomc/tbx3
  IPC1 = c(11,20,32), #Ascl1
  IPC2 = c(30,35,42), # Neurog2
  RGC = c(12,24,36,38,39),# cRGC
  AS = c(2,34),
  TC = c(4), #RGC/TC
  Microglia = 41, #Microglia/Fib
  OPC = 19,
  low_quality1 = 21,
  low_quality2 = 37,
  Ependymal = 40,
  OD = 43
)

obj$rpca_cell_type <- NA

# 根据映射字典为每个cluster添加cell_types列
for (i in seq_along(rpca_cell_type_mapping)) {
  clusters <- rpca_cell_type_mapping[[i]]
  cell_type <- names(rpca_cell_type_mapping)[i]
  
  # 只为 rpca_cell_type 为空的情况下进行赋值
  mask <- obj$rpca.clust %in% clusters & is.na(obj$rpca_cell_type)
  
  # 如果有满足条件的行，则进行赋值
  if (any(mask)) {
    obj$rpca_cell_type[mask] <- cell_type
  }
}



rpca_cell_type_dotplot<-DotPlot(obj, group.by = "rpca_cell_type", features = unique(unlist(marker_genes_list))) + RotatedAxis() +
  ggtitle("RPCA Celltype Markers") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

ggsave("rpca_cell_type_dotplot.png", rpca_cell_type_dotplot, width = 40, height = 8, dpi = 300)

VlnPlot(obj,features = c("nCount_RNA","nFeature_RNA","percent.mt"),group.by = "rpca_cell_type",pt.size = 0,log = T)

umap_split_by_batch<-DimPlot(
  obj,
  reduction = "rpca.umap",
  group.by = "rpca_cell_type",
  raster = F,
  split.by = "Batch",ncol = 4
) + ggtitle("integrated by RPCA")

ggsave("umap_split_by_batch.png", umap_split_by_batch, width = 20, height = 20, dpi = 300)


Neuron<-subset(obj,subset = rpca_cell_type == "GABA"| rpca_cell_type == "GLU"|rpca_cell_type == "Neu"|rpca_cell_type == "Neu_Pomc")

saveRDS(Neuron,"./Neuron.RDS")

DimPlot(
  obj,
  reduction = "rpca.umap",
  group.by = "rpca_cell_type",
  raster = F,
  label = T
) + ggtitle("RPCA celltype")

# Get unique cell types
cell_types <- unique(obj$rpca_cell_type)

# Create a list to store the plots
plot_list <- list()

# Generate a plot for each cell type
for (cell_type in cell_types) {
  plot <- DimPlot(
    obj,
    reduction = "rpca.umap",
    group.by = "rpca_cell_type",
    raster = FALSE,
    label = TRUE,
    cells.highlight = WhichCells(obj, expression = rpca_cell_type == cell_type),
    cols.highlight = "red",
    cols = "gray"
  ) + 
    ggtitle(paste("RPCA celltype -", cell_type)) +
    NoLegend()
  
  plot_list[[cell_type]] <- plot
}

# Combine all plots
combined_plot <- wrap_plots(plot_list, ncol = 4)  # Adjust ncol as needed

# Add a main title to the combined plot
final_plot <- combined_plot + 
  plot_annotation(title = "RPCA Celltypes Highlighted", 
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Display the combined plot
# print(final_plot)

# Save the combined plot
ggsave("RPCA_celltypes_celltype_combined.png", final_plot, width = 24, height = 24, limitsize = FALSE)

filtered_obj<-subset(obj,subset = rpca_cell_type != "low_quality1" & rpca_cell_type != "low_quality2")

# Get unique cell types
filtered_cell_types <- unique(filtered_obj$rpca_cell_type)

# Create a list to store the plots
filtered_plot_list <- list()

# Generate a plot for each cell type
for (cell_type in filtered_cell_types) {
  plot <- DimPlot(
    filtered_obj,
    reduction = "rpca.umap",
    group.by = "rpca_cell_type",
    raster = FALSE,
    label = TRUE,
    cells.highlight = WhichCells(filtered_obj, expression = rpca_cell_type == cell_type),
    cols.highlight = "red",
    cols = "gray"
  ) + 
    ggtitle(paste("RPCA celltype -", cell_type)) +
    NoLegend()
  
  filtered_plot_list[[cell_type]] <- plot
}

filtered_combined_plot<-NULL

# Combine all plots
filtered_combined_plot <- wrap_plots(filtered_plot_list, ncol = 4)  # Adjust ncol as needed

# Add a main title to the combined plot
filtered_final_plot <- filtered_combined_plot + 
  plot_annotation(title = "RPCA Celltypes", 
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Save the combined plot
ggsave("RPCA_filtered_celltypes.png", filtered_final_plot, width = 24, height = 18, limitsize = FALSE)

# 分时间点看umap图

# 提取发育时间点并创建新的列
filtered_obj$stage <- sub("_.*", "", obj$Batch)

filtered_obj@meta.data$percent.Malat1<-NULL
filtered_obj@meta.data$log10_nCount_RNA<-NULL
filtered_obj@meta.data$log10_nFeature_RNA<-NULL
filtered_obj@meta.data$log10_percent.hb<-NULL
filtered_obj@meta.data$log10_percent.mt<-NULL
filtered_obj@meta.data$log10_percent.ribo<-NULL
filtered_obj@meta.data$log10_percent.Malat1<-NULL 
filtered_obj@meta.data$log2_percent.mt<-NULL 

# save
saveRDS(filtered_obj,"./filtered_obj.rds")

filtered_obj<-readRDS("./filtered_jointpca_annotated_99252_cells.rds")

# Get unique cell types
stage_list <- unique(filtered_obj$stage)

# Create a list to store the plots
stage_plot_list <- list()

# Generate a plot for each cell type
for (dev_stage in stage_list) {
  plot <- DimPlot(
    filtered_obj,
    reduction = "rpca.umap",
    group.by = "rpca_cell_type",
    raster = FALSE,
    label = T,
    label.size = 5,
    cells.highlight = WhichCells(filtered_obj, expression = stage == dev_stage),
    cols.highlight = "red",
    cols = "gray"
  ) + 
    ggtitle(paste("RPCA stage -", dev_stage)) +
    NoLegend()
  
  stage_plot_list[[dev_stage]] <- plot
}

stage_combined_plot<-NULL

# Combine all plots
stage_combined_plot <- wrap_plots(stage_plot_list, ncol = 4)  # Adjust ncol as needed

# Add a main title to the combined plot
stage_final_plot <- stage_combined_plot + 
  plot_annotation(title = "RPCA stages", 
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Save the combined plot
ggsave("RPCA_stages.png", stage_final_plot, width = 24, height = 12, limitsize = FALSE)

# DimPlot(
#   obj,
#   reduction = "rpca.umap",
#   group.by = "rpca_cell_type",
#   raster = F,
#   split.by = "rpca_cell_type",ncol = 4
# ) + ggtitle("integrated by RPCA")

library(Nebulosa)
plot_density(obj,marker_genes_list$RGC,reduction = "rpca.umap")
plot_density(obj,marker_genes_list$IPC,reduction = "rpca.umap")
plot_density(obj,c("Olig1","Pdgfra","Mbp","Mag"),reduction = "rpca.umap")
plot_density(obj,c(marker_genes_list$GABA,marker_genes_list$GLU),reduction = "rpca.umap")
plot_density(obj,c(marker_genes_list$Fibroblasts),reduction = "rpca.umap")
plot_density(obj,c(marker_genes_list$Microglia),reduction = "rpca.umap")


# Harmony


# SCVI-------------------------------------------------------

# UMAP降维查看批次效应
obj <- RunUMAP(
  obj,
  dims = 1:15,
  n.neighbors = 100,
  spread = 3,
  min.dist = 2,#min_dist must be less than or equal to spread
  reduction = "integrated.scvi",
  reduction.name = "scvi.umap"
)

DimPlot(obj, group.by = "Batch", reduction = "scvi.umap") + ggtitle("integrated by SCVI")


# 聚类
obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.scvi")
obj <- FindClusters(
  obj,
  cluster.name = "scvi.clust",
  method = "igraph",
  algorithm = 4,
  resolution = 0.4,
  verbose = T
)

# 可视化聚类结果
DimPlot(obj, reduction = "scvi.umap", group.by = "scvi.clust")

DotPlot(obj, group.by = "scvi.clust", features = marker_genes_list) + RotatedAxis() +
  ggtitle("integrated by SCVI") +
  theme(plot.title = element_text(hjust = 0.5))

# 创建一个细胞类型到cluster的映射字典 25不确定
rpca_cell_type_mapping <- list(
  GABA = c(2, 3, 6, 9, 10, 19, 20, 21, 24),
  GLU = c(1, 5, 11, 12, 15, 16, 17, 23),
  IPC_Ascl1 = c(8, 13),
  IPC_Neurog2 = c(22),
  TC_RGC = c(14),
  AS_RGC = c(4),
  RGC = c(7),
  Microglia_Fib = 27,
  OD_OPC = 18,
  EC = 26,
  low_quality = 25
)

obj$rpca_cell_type <- NA

# 根据映射字典为每个cluster添加cell_types列
for (i in seq_along(rpca_cell_type_mapping)) {
  clusters <- rpca_cell_type_mapping[[i]]
  cell_type <- names(rpca_cell_type_mapping)[i]
  obj$rpca_cell_type[obj$rpca.clust %in% clusters] <- cell_type
}

obj <- RunUMAP(
  obj,
  dims = 1:30,
  n.neighbors = 35,
  spread = 0.25,
  min.dist = 0.25,#min_dist must be less than or equal to spread
  reduction = "integrated.rpca",
  reduction.name = "rpca.umap",
  umap.method = "umap-learn",
  metric = "euclidean"
)

DimPlot(
  obj,
  reduction = "rpca.umap",
  group.by = "rpca_cell_type",
  raster = F,
  label = T
) + ggtitle("integrated by RPCA")

DimPlot(
  obj,
  reduction = "rpca.umap",
  group.by = "rpca_cell_type",
  raster = F,
  split.by = "rpca_cell_type",ncol = 4
) + ggtitle("integrated by RPCA")


# COSG
remotes::install_github(repo = 'genecell/COSGR')
library(COSG)

obj<-JoinLayers(obj)

COSG_markers <- cosg(
  obj,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1,
  remove_lowly_expressed = TRUE,
  expressed_pct = 0.25,
  n_genes_user = 50
)

DotPlot(obj, unique(unlist(COSG_markers$names)))+RotatedAxis()

# 输出结果
# 将向量转换为数据框
df_COSG_markers <- data.frame(Marker = COSG_markers$names,
                              Value = COSG_markers$scores,
                              stringsAsFactors = FALSE)

# 重新排列列
n_columns <- ncol(df_COSG_markers) / 2
column_order <- as.vector(rbind(1:n_columns, (n_columns+1):(2*n_columns)))

df_COSG_markers <- df_COSG_markers[, column_order]

# 保存为 CSV 文件
write.csv(df_COSG_markers, file = "COSG_markers.csv", row.names = FALSE)


