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

# scvi
obj <- IntegrateLayers(
  object = obj,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/opt/miniforge3/envs/scvi/",
  verbose = FALSE
)

gc()

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
