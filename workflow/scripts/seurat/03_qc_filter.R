library(Seurat)
library(ggplot2)

# 从命令行参数获取输入、输出和配置路径
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
qc_params <- snakemake@params[["qc"]]

# 打印当前工作目录和传递的参数
print(paste("当前工作目录:", getwd()))
print(paste("输入文件路径:", input_file))
print(paste("输出文件路径:", output_file))
print(paste("QC 参数:", toString(qc_params)))

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)

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

