library(Seurat)
# BiocManager::install("gprofiler2")
library(gprofiler2)
# 从 Seurat 获取细胞周期基因
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

# 转换基因名
mmus_s <- gorth(s_genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m <- gorth(g2m_genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# 保存转换后的基因名到文件
write.table(mmus_s, file="mmus_s_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(mmus_g2m, file="mmus_g2m_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

