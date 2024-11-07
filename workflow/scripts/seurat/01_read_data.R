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

