library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(scCustomize)

# 从命令行参数获取输入、输出和配置路径
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
figures_dir <- snakemake@params[["figures_dir"]]
n_var_features <- snakemake@params[["n_var_features"]]
n_pcs <- snakemake@params[["n_pcs"]]
percent_doublets <- snakemake@params[["percent_doublets"]]
regress_features <- snakemake@params[["regress_features"]]

# 打印当前工作目录和传递的参数
print(paste("当前工作目录:", getwd()))
print(paste("输入文件路径:", input_file))
print(paste("输出文件路径:", output_file))
print(paste("图表输出目录:", figures_dir))
print(paste("高变基因数量:", n_var_features))
print(paste("PCA 组件数:", n_pcs))
print(paste("双重体百分比:", percent_doublets))
print(paste("回归特征:", paste(regress_features, collapse = ", ")))

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)

# obj <- JoinLayers(obj)

# 细胞周期基因读入
genes.df <- read.csv("workflow/scripts/seurat/genes.csv", stringsAsFactors = FALSE)

# 计算细胞周期评分
obj <- CellCycleScoring(
  obj,
  s.features = genes.df$s.genes,
  g2m.features = genes.df$g2m.genes,
  set.ident = TRUE
)

# 分割数据按 Batch
seurat_list <- SplitObject(obj, split.by = "Batch")
rm(obj)
gc()

# 预处理数据：NormalizeData, FindVariableFeatures, ScaleData, RunPCA
pipeline <- function(obj, n_var_features, PCA_npcs, regress_features) {
  obj <- obj %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = n_var_features) %>%
    ScaleData(regress.vars = regress_features) %>% # 使用回归特征
    RunPCA(npcs = PCA_npcs, verbose = FALSE) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters()
  return(obj)
}

# 运行预处理管道
seurat_list <- lapply(seurat_list, pipeline, n_var_features = n_var_features, PCA_npcs = n_pcs, regress_features = regress_features)

# 运行 DoubletFinder
Find_doublet <- function(data, PCA_npcs, Percent_doublets) {
  sweep.res.list <- paramSweep(data, PCs = 1:PCA_npcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(Percent_doublets * ncol(data))
  p <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK)) # pK Selection
  
  # Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi <- round(Percent_doublets * nrow(data@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # 运行 DoubletFinder
  data <- doubletFinder(
    data,
    PCs = 1:PCA_npcs,
    pN = 0.25,
    pK = p,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  # 保存双重体信息到 meta.data
  colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
  return(data)
}

# 对每个 Batch 运行 DoubletFinder
seurat_list <- lapply(seurat_list, function(x) {
  Find_doublet(x, PCA_npcs = n_pcs, Percent_doublets = percent_doublets)
})

# 合并所有结果
all_merged_postQC_add_doublet_info <- Merge_Seurat_List(seurat_list)

# 保存合并后的结果
saveRDS(all_merged_postQC_add_doublet_info, output_file)

print("数据处理、PCA和DoubletFinder成功完成")

