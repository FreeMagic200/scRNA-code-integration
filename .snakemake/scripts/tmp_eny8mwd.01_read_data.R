
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
    input = list('data/raw/seurat_input.rds', "rds" = 'data/raw/seurat_input.rds'),
    output = list('data/processed/seurat/01_read_data.rds', "rds" = 'data/processed/seurat/01_read_data.rds'),
    params = list('data/raw/seurat_input.rds', 'data/processed/seurat/01_read_data.rds', "input_file" = 'data/raw/seurat_input.rds', "output_file" = 'data/processed/seurat/01_read_data.rds'),
    wildcards = list(),
    threads = 1,
    log = list('logs/seurat/01_read_data.log'),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("input" = list("seurat" = 'data/raw/seurat_input.rds', "scanpy" = 'data/raw/scanpy_input.h5ad'), "output_dir" = list("seurat" = 'results/seurat', "scanpy" = 'results/scanpy'), "processed_dir" = list("seurat" = 'data/processed/seurat', "scanpy" = 'data/processed/scanpy'), "figures_dir" = list("seurat" = 'figures/seurat', "scanpy" = 'figures/scanpy'), "qc" = list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), "n_var_features" = 2500, "regress_features" = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score'), "n_pcs" = 30, "resolutions" = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3), "umap" = list("dims" = 30, "n_neighbors" = 30, "min_dist" = 0.4, "spread" = 1, "method" = 'umap-learn', "metric" = 'minkowski'), "threads" = 1, "random_seed" = 42, "logging" = list("level" = 'INFO', "output" = 'logs/analysis.log')),
    rule = 'seurat_read_data',
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
# 加载必要的库
library(Seurat)

cat("当前工作目录: ", getwd(), "\n")

# 获取命令行参数
input_file<-snakemake@params[["input_file"]]
output_file<-snakemake@params[["output_file"]]

# 打印所有传递的命令行参数
cat("传递的参数:\n")
print(input_file)

# 定义输入文件路径
cat("正在从以下路径读取 Seurat 对象:", input_file, "\n")

tryCatch({
  # 从 RDS 文件读取 Seurat 对象
  obj <- readRDS(input_file)
  cat("成功从以下路径读取 Seurat 对象:", input_file, "\n")
  
  # 更新 Seurat 对象
  obj <- UpdateSeuratObject(obj)
  cat("Seurat 对象成功更新\n")
  
  # 检查对象完整性
  if (!"Seurat" %in% class(obj)) {
    stop("对象不是 'Seurat' 类")
  }
  
  # 其他完整性检查示例
  if (is.null(obj@assays) || length(obj@assays) == 0) {
    stop("Seurat 对象中没有 assays")
  }
  if (is.null(obj@meta.data) || nrow(obj@meta.data) == 0) {
    stop("Seurat 对象中没有 metadata")
  }
  
  cat("Seurat 对象通过完整性检查\n")
  
}, error = function(e) {
  # 捕获错误并记录错误信息
  cat("发生错误:", e$message, "\n")
  quit(status = 1)
})

saveRDS(obj,paste0(output_file))

cat("脚本执行成功\n")

