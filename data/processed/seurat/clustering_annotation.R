library(Seurat)
library(ggplot2)
library(stringr)
library(future)
library(patchwork)
library(clustree)
library(dplyr)
library(scCustomize)
set.seed(42)

obj<-readRDS("./umap_tunned_103530_cells.rds")
plan("sequential", workers = 1)
options(future.globals.maxSize = Inf, future.seed = TRUE,future.seed = T)
plan()


# jointPCA-------------------------------------------------------
DimPlot(obj, group.by = "Batch", reduction = "jointpca.umap") + ggtitle("integrated by jointpca")

jointpca_clust_umap<-DimPlot(
  unobj,
  reduction = "jointpca.umap",
  group.by = "jointpca.clust",
  raster = F,
  label = T
) + ggtitle("UMAP of jointpca clusters")

ggsave("jointpca_clust_umap.png",plot = jointpca_clust_umap,width = 18.15,height = 14.34,units = "cm")


# 聚类
# 细分聚类分出OD
obj <- FindNeighbors(obj, dims = 1:12, reduction = "integrated.jointpca",k.param = 5,annoy.metric = "euclidean")
obj <- FindClusters(
  obj,
  method = "igraph",
  algorithm = 4,
  resolution = seq(0,2,0.1),
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
ggsave("k5_clustree_indent_1.png", plot = k5_tree, width = 12, height = 14, dpi = 200)

DimPlot(obj, reduction = "jointpca.umap", group.by = "RNA_snn_res.1.1",raster = F,label = T)

# # 转换为数值型，并重新设置因子水平
obj$RNA_snn_res.1.1 <- factor(obj$RNA_snn_res.1.1, levels = as.character(sort(as.numeric(levels(obj$RNA_snn_res.1.1)))))

# Step 1: 展平列表为一个向量
flattened_genes <- unlist(marker_genes_list)

# Step 2: 去除重复项
unique_genes <- unique(flattened_genes)

# Step 3: 根据原始列表结构重建新的列表
reconstructed_list <- list()

for (gene in unique_genes) {
  # 找到第一个包含该基因的组名
  group_name <- names(marker_genes_list)[sapply(marker_genes_list, function(x) gene %in% x)][1]
  
  # 如果该组已经存在于新列表中，则添加该基因
  if (group_name %in% names(reconstructed_list)) {
    reconstructed_list[[group_name]] <- c(reconstructed_list[[group_name]], gene)
  } else {
    # 如果该组不存在，则创建新的组并添加该基因
    reconstructed_list[[group_name]] <- gene
  }
}

# 输出结果
reconstructed_list

# 创建点图
p_1.1 <- DotPlot(obj, group.by = "RNA_snn_res.1.1", features = reconstructed_list) + 
  RotatedAxis() +
  ggtitle("integrated by jointPCA") +
  theme(plot.title = element_text(hjust = 0.5))

# 保存图形
ggsave("res1.1_jointPCA_clust.png", plot = p_1.1, width = 35, height = 12, dpi = 300)

obj$jointpca.clust <-obj$RNA_snn_res.1.1

# 转换为数值型，并重新设置因子水平
obj$jointpca.clust <- factor(obj$jointpca.clust, levels = as.character(sort(as.numeric(levels(obj$jointpca.clust)))))

FeaturePlot(obj,c("Crym","Col25a1","Rax"),order = F,raster = F,max.cutoff = 1.15,pt.size = 0.1)

# 获取基因列表
genes_to_plot <- c(unlist(reconstructed_list), "Col25a1")

# 计算基因数量
num_genes <- length(genes_to_plot)

# 动态调整 ncol
ncol <- min(5, ceiling(sqrt(num_genes)))  # 最多5列，但会根据基因数量调整

# 创建 FeaturePlot
p <- FeaturePlot(obj, 
                 features = genes_to_plot, 
                 order = FALSE, 
                 raster = FALSE, 
                 max.cutoff = 1.18,
                 pt.size = 0.1, 
                 ncol = ncol,
                 reduction = "jointpca.umap")

p1 <- FeaturePlot(obj, 
                  features = genes_to_plot, 
                  order = T, 
                  raster = FALSE, 
                  max.cutoff = 1.18,
                  pt.size = 0.1, 
                  ncol = ncol,
                  reduction = "jointpca.umap")

# 计算合理的图片大小
width <- ncol * 6  # 每列5英寸宽
height <- ceiling(num_genes / ncol) * 5  # 每行5英寸高

# 保存图片
ggsave("FeaturePlot_genes.png", plot = p, width = width, height = height, dpi = 300, limitsize = FALSE)

ggsave("FeaturePlot_ordered_genes.png", plot = p1, width = width, height = height, dpi = 300, limitsize = FALSE)

VlnPlot(obj,c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb","Malat1"),group.by = "RNA_snn_res.1.1",pt.size = 0,ncol = 2,log = T)

jointpca_clust_dotplot<-DotPlot(obj, group.by = "RNA_snn_res.1.1", features = reconstructed_list) + RotatedAxis() +
  ggtitle("JointPCA clust Markers") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

ggsave("jointpca_clust_dotplot.png", jointpca_clust_dotplot, width = 34, height = 15, dpi = 300)

DimPlot(obj, reduction = "jointpca.umap", group.by = "RNA_snn_res.1.1",raster = F,label = T)

# 寻找所有聚类的标记基因
obj <- JoinLayers(obj)
scRNA.markers <-
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.6,
    slot = "data"
  )

write.csv(scRNA.markers, file = "./scRNA_all_markers.csv", row.names = FALSE)

all.markers<-read.csv("./scRNA_all_markers.csv")

# 先筛选 avg_log2FC > 0，再筛选 p_val < 0.05
all.markers <- all.markers %>%
  dplyr::select(gene, everything()) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::filter(p_val < 0.05)

# 将avg_log2FC排名前50的基因筛选出来
top50 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

write.csv(top50, file = "./scRNA_top50_markers.csv", row.names = FALSE)

top50<-read.csv("./scRNA_top50_markers.csv")

# 创建一个细胞类型到cluster的映射字典 25不确定
jointpca_cell_type_mapping <- list(
  GABA = c(1,2,7,8,9,13,21,23,25),
  GLU = c(3,4,5,6,14,16,18,24,33,35), #24很少表达Slc17a6/Slc32a1 14/35都表达
  IPC1 = c(12,20), #Ascl1
  IPC2 = c(28,29), # Neurog2
  RGC = c(11,15,26,27,31),# cRGC
  AS = c(17,22),
  TC = c(30), #RGC/TC
  Mig = 34, #Microglia/Fib
  OPC = 19,
  low_quality1 = 10,
  low_quality2 = 32,
  Ependymal = 37,
  OD = 36
)

obj$jointpca_cell_type <- NA

# 根据映射字典为每个cluster添加cell_types列
for (i in seq_along(jointpca_cell_type_mapping)) {
  clusters <- jointpca_cell_type_mapping[[i]]
  cell_type <- names(jointpca_cell_type_mapping)[i]
  
  # 只为 jointpca_cell_type 为空的情况下进行赋值
  mask <- obj$jointpca.clust %in% clusters & is.na(obj$jointpca_cell_type)
  
  # 如果有满足条件的行，则进行赋值
  if (any(mask)) {
    obj$jointpca_cell_type[mask] <- cell_type
  }
}

jointpca_cell_type_dotplot<-DotPlot(obj, group.by = "jointpca_cell_type", features = unique(unlist(marker_genes_list))) + RotatedAxis() +
  ggtitle("jointpca Celltype Markers") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

ggsave("jointpca_cell_type_dotplot.png", jointpca_cell_type_dotplot, width = 40, height = 8, dpi = 300)

cell_type_QC_vln_plot<-VlnPlot(obj,c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.hb","Malat1"),group.by = "jointpca_cell_type",pt.size = 0,ncol = 2,log = F)

ggsave("cell_type_QC_vln_plot.png", width = 28, height = 18, units = "cm")

umap_split_by_batch<-DimPlot(
  obj,
  reduction = "jointpca.umap",
  group.by = "jointpca_cell_type",
  raster = F,
  split.by = "Batch",ncol = 4
) + ggtitle("integrated by jointpca")

ggsave("umap_split_by_batch.png", umap_split_by_batch, width = 20, height = 20, dpi = 300)

saveRDS(obj,"./jointpca_annotated_103530_cells.rds")

obj<-readRDS("./jointpca_annotated_103530_cells.rds")

jointpca_unfiltered_cell_type_umap<-DimPlot(
  obj,
  reduction = "jointpca.umap",
  group.by = "jointpca_cell_type",
  raster = F,
  label = T
) + ggtitle("jointpca celltype")

ggsave("jointpca_unfiltered_cell_type_umap.png",plot = jointpca_unfiltered_cell_type_umap,width = 18.15,height = 14.34,units = "cm")

obj<-subset(obj,subset = jointpca_cell_type != "low_quality1" & jointpca_cell_type != "low_quality2")

saveRDS(obj,"./filtered_jointpca_annotated_99252_cells.rds")

obj<-readRDS("./filtered_jointpca_annotated_99252_cells.rds")

# 寻找所有聚类的标记基因
obj<-JoinLayers(obj)
filtered_scRNA.markers <-
  FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.6,
    slot = "data"
  )

write.csv(filtered_scRNA.markers, file = "./filtered_scRNA.markers.csv", row.names = FALSE)

filtered_scRNA.markers<-read.csv("./filtered_scRNA.markers.csv")

Neuron<-subset(obj,subset = jointpca_cell_type == "GABA"| jointpca_cell_type == "GLU")

saveRDS(Neuron,"./Neuron.RDS")

jointpca_cell_type_umap<-DimPlot(
  obj,
  reduction = "jointpca.umap",
  group.by = "jointpca_cell_type",
  raster = F,
  label = T
) + ggtitle("jointpca celltype")

ggsave("jointpca_cell_type_umap.png",width = 18.15,height = 14.34,units = "cm")

# Get unique cell types
cell_types <- unique(obj$jointpca_cell_type)

# Create a list to store the plots
plot_list <- list()

# Generate a plot for each cell type
for (cell_type in cell_types) {
  plot <- DimPlot(
    obj,
    reduction = "jointpca.umap",
    group.by = "jointpca_cell_type",
    raster = FALSE,
    label = TRUE,
    cells.highlight = WhichCells(obj, expression = jointpca_cell_type == cell_type),
    cols.highlight = "red",
    cols = "gray"
  ) + 
    ggtitle(paste("jointpca celltype -", cell_type)) +
    NoLegend()
  
  plot_list[[cell_type]] <- plot
}

# Combine all plots
combined_plot <- wrap_plots(plot_list, ncol = 4)  # Adjust ncol as needed

# Add a main title to the combined plot
final_plot <- combined_plot + 
  plot_annotation(title = "jointpca Celltypes Highlighted", 
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Display the combined plot
# print(final_plot)

# Save the combined plot
ggsave("jointpca_celltypes_celltype_combined.png", final_plot, width = 30, height = 24, limitsize = FALSE)

# Get unique cell types
filtered_cell_types <- unique(obj$jointpca_cell_type)

# Create a list to store the plots
filtered_plot_list <- list()

# 分时间点看umap图

# 提取发育时间点并创建新的列
obj$stage <- sub("_.*", "", obj$Batch)

obj@meta.data$percent.Malat1<-NULL
obj@meta.data$log10_nCount_RNA<-NULL
obj@meta.data$log10_nFeature_RNA<-NULL
obj@meta.data$log10_percent.hb<-NULL
obj@meta.data$log10_percent.mt<-NULL
obj@meta.data$log10_percent.ribo<-NULL
obj@meta.data$log10_percent.Malat1<-NULL 
obj@meta.data$log2_percent.mt<-NULL 

# save
saveRDS(obj,"./filtered_jointpca_annotated_99252_cells.rds")

obj<-readRDS("./filtered_jointpca_annotated_99252_cells.rds")

# Get unique cell types
stage_list <- unique(obj$stage)

# Create a list to store the plots
stage_plot_list <- list()

# Generate a plot for each cell type
for (dev_stage in stage_list) {
  plot <- DimPlot(
    obj,
    reduction = "jointpca.umap",
    group.by = "jointpca_cell_type",
    raster = FALSE,
    label = T,
    label.size = 5,
    cells.highlight = WhichCells(obj, expression = stage == dev_stage),
    cols.highlight = "red",
    cols = "gray"
  ) + 
    ggtitle(paste("jointpca stage -", dev_stage)) +
    NoLegend()
  
  stage_plot_list[[dev_stage]] <- plot
}

stage_combined_plot<-NULL

# Combine all plots
stage_combined_plot <- wrap_plots(stage_plot_list, ncol = 4)  # Adjust ncol as needed

# Add a main title to the combined plot
stage_final_plot <- stage_combined_plot + 
  plot_annotation(title = "jointpca stages", 
                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Save the combined plot
ggsave("jointpca_stages.png", stage_final_plot, width = 24, height = 12, limitsize = FALSE)

# DimPlot(
#   obj,
#   reduction = "jointpca.umap",
#   group.by = "jointpca_cell_type",
#   raster = F,
#   split.by = "jointpca_cell_type",ncol = 4
# ) + ggtitle("integrated by jointpca")

library(Nebulosa)
plot_density(obj,marker_genes_list$RGC,reduction = "jointpca.umap")
plot_density(obj,marker_genes_list$IPC,reduction = "jointpca.umap")
plot_density(obj,c("Olig1","Pdgfra","Mbp","Mag"),reduction = "jointpca.umap")
plot_density(obj,c(marker_genes_list$GABA,marker_genes_list$GLU),reduction = "jointpca.umap")
plot_density(obj,c(marker_genes_list$Fibroblasts),reduction = "jointpca.umap")
plot_density(obj,c(marker_genes_list$Microglia),reduction = "jointpca.umap")






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



