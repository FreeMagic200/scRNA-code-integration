
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
    input = list('data/processed/seurat/02_normalized_data.rds', "rds" = 'data/processed/seurat/02_normalized_data.rds'),
    output = list('data/processed/seurat/03_filtered_data.rds', "rds" = 'data/processed/seurat/03_filtered_data.rds'),
    params = list(list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), 'figures/seurat', "qc" = list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), "figures_dir" = 'figures/seurat'),
    wildcards = list(),
    threads = 1,
    log = list('logs/seurat/03_qc_filter.log'),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("input" = list("seurat" = 'data/raw/seurat_input.rds', "scanpy" = 'data/raw/scanpy_input.h5ad'), "output_dir" = list("seurat" = 'results/seurat', "scanpy" = 'results/scanpy'), "processed_dir" = list("seurat" = 'data/processed/seurat', "scanpy" = 'data/processed/scanpy'), "figures_dir" = list("seurat" = 'figures/seurat', "scanpy" = 'figures/scanpy'), "qc" = list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), "n_var_features" = 2500, "regress_features" = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score'), "n_pcs" = 30, "resolutions" = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3), "umap" = list("dims" = 30, "n_neighbors" = 30, "min_dist" = 0.4, "spread" = 1, "method" = 'umap-learn', "metric" = 'minkowski'), "threads" = 1, "random_seed" = 42, "logging" = list("level" = 'INFO', "output" = 'logs/analysis.log')),
    rule = 'seurat_qc_filter',
    bench_iteration = as.numeric(NA),
    scriptdir = '/mnt/data/projects/scRNA-0722/workflow/scripts/seurat',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
library(Seurat)
library(ggplot2)

# 从命令行参数获取输入、输出和配置路径
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
qc_params <- snakemake@params[["qc"]]
figures_dir <- snakemake@params[["figures_dir"]]

# 打印当前工作目录和传递的参数
print(paste("当前工作目录:", getwd()))
print(paste("输入文件路径:", input_file))
print(paste("输出文件路径:", output_file))
print(paste("QC 参数:", toString(qc_params)))

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)

# 计算 QC 指标
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(egf)|^(s1l)|^(p1)].+")

# 绘制 QC 图
pdf(file = file.path(figures_dir, "qc_plots.pdf"))
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 2)
dev.off()

# QC 过滤函数
QC_filter <- function(obj, qc_params) {
  obj <- subset(
    obj,
    subset = nFeature_RNA < qc_params$max_nFeature_RNA &
      nFeature_RNA > qc_params$min_nFeature_RNA &
      nCount_RNA < qc_params$max_nCount_RNA &
      nCount_RNA > qc_params$min_nCount_RNA &
      percent.mt < qc_params$max_percent_mt &
      percent.ribo < qc_params$max_percent_ribo &
      Malat1 > qc_params$min_num_Malat1 &
      percent.hb < qc_params$max_percent_hb
  )
  return(obj)
}

# 应用 QC 过滤
obj_filtered <- QC_filter(obj, qc_params)

# 保存过滤后的 Seurat 对象
saveRDS(obj_filtered, output_file)

# 打印成功信息
print("QC 处理和过滤成功")

