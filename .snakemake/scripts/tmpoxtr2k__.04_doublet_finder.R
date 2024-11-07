
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('data/processed/seurat/03_filtered_data.rds', "rds" = 'data/processed/seurat/03_filtered_data.rds'),
    output = list('data/processed/seurat/04_doublet_data.rds', "rds" = 'data/processed/seurat/04_doublet_data.rds'),
    params = list('figures/seurat', 2500, 30, 0.05, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score'), "figures_dir" = 'figures/seurat', "n_var_features" = 2500, "n_pcs" = 30, "percent_doublets" = 0.05, "regress_features" = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score')),
    wildcards = list(),
    threads = 1,
    log = list('logs/seurat/04_doublet_finder.log'),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("input" = list("seurat" = 'data/raw/seurat_input.rds', "scanpy" = 'data/raw/scanpy_input.h5ad'), "output_dir" = list("seurat" = 'results/seurat', "scanpy" = 'results/scanpy'), "processed_dir" = list("seurat" = 'data/processed/seurat', "scanpy" = 'data/processed/scanpy'), "figures_dir" = list("seurat" = 'figures/seurat', "scanpy" = 'figures/scanpy'), "qc" = list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), "n_var_features" = 2500, "regress_features" = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score'), "n_pcs" = 30, "resolutions" = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3), "umap" = list("dims" = 30, "n_neighbors" = 30, "min_dist" = 0.4, "spread" = 1, "method" = 'umap-learn', "metric" = 'minkowski'), "threads" = 1, "random_seed" = 42, "logging" = list("level" = 'INFO', "output" = 'logs/analysis.log')),
    rule = 'seurat_doublet_finder',
    bench_iteration = as.numeric(NA),
    scriptdir = '/mnt/data/projects/scRNA-0722-Reviewed/workflow/scripts/seurat',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
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

# 细胞周期基因评分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 转换为小写
s.genes <- paste0(substr(s.genes, 1, 1), tolower(substr(s.genes, 2, nchar(s.genes))))
g2m.genes <- paste0(substr(g2m.genes, 1, 1), tolower(substr(g2m.genes, 2, nchar(g2m.genes))))

# 计算细胞周期评分
obj <- CellCycleScoring(
  obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# 分割数据按 Batch
seurat_list <- SplitObject(obj, split.by = "Batch")
rm(obj)
gc()

# 预处理数据：NormalizeData, FindVariableFeatures, ScaleData, RunPCA
pipeline <- function(obj, n_var_features, PCA_npcs, regress_features) {
  obj <- obj %>%
    # NormalizeData()%>%
    FindVariableFeatures(selection.method = "vst", nfeatures = n_var_features) %>%
    ScaleData(regress.vars = regress_features) %>% # 使用回归特征
    RunPCA(npcs = PCA_npcs, verbose = FALSE)
  return(obj)
}

# 运行预处理管道
seurat_list <- lapply(seurat_list, pipeline, n_var_features = n_var_features, PCA_npcs = n_pcs, regress_features = regress_features)

# 运行 DoubletFinder
Find_doublet <- function(data, PCA_npcs, Percent_doublets) {
  sweep.res.list <- paramSweep(data, PCs = 1:PCA_npcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK)) # pK Selection
  
  # Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi <- round(Percent_doublets * ncol(data))
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

