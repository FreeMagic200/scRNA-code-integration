library(Seurat)

# 从命令行参数获取输入、输出和配置路径
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
n_var_features <- snakemake@params[["n_var_features"]]
n_pcs <- snakemake@params[["n_pcs"]]
regress_features <- snakemake@params[["regress_features"]]

# 打印当前工作目录和传递的参数
print(paste("当前工作目录:", getwd()))
print(paste("输入文件路径:", input_file))
print(paste("输出文件路径:", output_file))
print(paste("高变基因数量:", n_var_features))
print(paste("PCA 组件数:", n_pcs))
print(paste("回归特征:", paste(regress_features, collapse = ", ")))

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)
# 再次预处理数据
obj <- subset(obj,doublet_info == "Singlet")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_var_features)
obj <- ScaleData(obj, vars.to.regress = regress_features)
obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)

# 保存后处理结果
saveRDS(obj, output_file)

print("数据处理和PCA成功完成")

