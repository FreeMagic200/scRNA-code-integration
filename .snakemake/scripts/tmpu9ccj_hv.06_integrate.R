
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
    input = list('data/processed/seurat/05_postQC_data.rds', "rds" = 'data/processed/seurat/05_postQC_data.rds'),
    output = list('data/processed/seurat/06_integrated_data.rds', "rds" = 'data/processed/seurat/06_integrated_data.rds'),
    params = list(),
    wildcards = list(),
    threads = 1,
    log = list('logs/seurat/06_integrate.log'),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("input" = list("seurat" = 'data/raw/seurat_input.rds', "scanpy" = 'data/raw/scanpy_input.h5ad'), "output_dir" = list("seurat" = 'results/seurat', "scanpy" = 'results/scanpy'), "processed_dir" = list("seurat" = 'data/processed/seurat', "scanpy" = 'data/processed/scanpy'), "figures_dir" = list("seurat" = 'figures/seurat', "scanpy" = 'figures/scanpy'), "qc" = list("max_nFeature_RNA" = 8000, "min_nFeature_RNA" = 500, "max_nCount_RNA" = 30000, "min_nCount_RNA" = 1000, "max_percent_mt" = 10, "max_percent_ribo" = 20, "min_num_Malat1" = 3, "max_percent_hb" = 1, "percent_doublets" = 0.05), "n_var_features" = 2500, "regress_features" = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'G2M_score', 'S_score'), "n_pcs" = 30, "resolutions" = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3), "umap" = list("dims" = 30, "n_neighbors" = 30, "min_dist" = 0.4, "spread" = 1, "method" = 'umap-learn', "metric" = 'minkowski'), "threads" = 1, "random_seed" = 42, "logging" = list("level" = 'INFO', "output" = 'logs/analysis.log')),
    rule = 'seurat_integrate',
    bench_iteration = as.numeric(NA),
    scriptdir = '/mnt/data/projects/scRNA-0722-ReQC/workflow/scripts/seurat',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
set.seed(42)
library(Seurat)
library(dplyr)
library(stringr)
library(SeuratWrappers)

# 从命令行参数获取输入、输出和配置路径
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# 打印当前工作目录和传递的参数
print(paste("当前工作目录:", getwd()))
print(paste("输入文件路径:", input_file))
print(paste("输出文件路径:", output_file))

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)

gc()

options(future.globals.maxSize = Inf)

# Split the Seurat object by Batch
# obj[["RNA"]] <- SplitObject(obj[["RNA"]], split.by = "Batch")

# Perform JointPCA integration
obj <- IntegrateLayers(
  object = obj,
  method = JointPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.jointpca",
  verbose = FALSE
)

gc()

# mnn
obj <- IntegrateLayers(
  object = obj,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

gc()

# Perform harmony integration
obj <- IntegrateLayers(
  object = obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  verbose = FALSE
)

gc()

# scvi
obj <- IntegrateLayers(
  object = obj,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/opt/miniforge3/envs/scvi/bin/python",
  verbose = FALSE
)

gc()

# Perform CCA integration
obj <- IntegrateLayers(
  object = obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = TRUE
)

gc()

# Perform RPCA integration
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = TRUE
)

gc()

# 保存整合后的结果
saveRDS(obj, output_file)

print("数据整合成功完成")
