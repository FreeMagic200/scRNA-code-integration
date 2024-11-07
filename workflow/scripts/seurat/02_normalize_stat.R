library(Seurat)

cat("当前工作目录: ", getwd(), "\n")

# 获取命令行参数
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
figures_dir <- snakemake@params[["figures_dir"]]

# 打印所有传递的命令行参数
cat("传递的参数:\n")
cat("输入文件:", input_file, "\n")
cat("输出文件:", output_file, "\n")

# 从 RDS 文件读取 Seurat 对象
obj <- readRDS(input_file)
cat("成功从以下路径读取 Seurat 对象:", input_file, "\n")

# 计算 QC 指标
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(egf)|^(s1l)|^(p1)].+")

# 绘制 QC 图
pdf(file = file.path(figures_dir, "qc_plots.pdf"))
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), ncol = 2,pt.size = 0)
dev.off()

# 正规化数据
obj <- NormalizeData(obj)
cat("数据规范化成功\n")

saveRDS(obj, output_file)
cat("成功保存 Seurat 对象到: ", output_file, "\n")
cat("脚本执行成功\n")

