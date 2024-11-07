# 细胞周期基因评分
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 转换为小写
s.genes <- paste0(substr(s.genes, 1, 1), tolower(substr(s.genes, 2, nchar(s.genes))))
g2m.genes <- paste0(substr(g2m.genes, 1, 1), tolower(substr(g2m.genes, 2, nchar(g2m.genes))))

# 将s.genes和g2m.genes合并为数据框
genes.df <- data.frame(s.genes = s.genes, g2m.genes = g2m.genes)

# 将s.genes和g2m.genes分别作为数据框的两列
genes.df <- data.frame(
  s.genes = c(s.genes, rep(NA, max(length(s.genes), length(g2m.genes)) - length(s.genes))),
  g2m.genes = c(g2m.genes, rep(NA, max(length(s.genes), length(g2m.genes)) - length(g2m.genes)))
)

# 保存为CSV文件
write.csv(genes.df, file = "genes.csv", row.names = FALSE)

# 已经核查同源性