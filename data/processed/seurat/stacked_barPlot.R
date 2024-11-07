library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
set.seed(42)

# obj <- readRDS("./filtered_obj.rds")

# 生成发育时间点/细胞类型的计数表
stage_celltype_count <- table(obj$stage, obj$jointpca_cell_type)

# 计算每个细胞类型的总细胞数
celltype_totals <- colSums(stage_celltype_count)

# 计算每个发育时间点在每种细胞类型中的比例
stage_celltype_proportion <- as.data.frame(t(t(stage_celltype_count) / celltype_totals))

# 将数据转换为长格式
stage_celltype_proportion_long <- melt(stage_celltype_proportion)

stage_celltype_proportion_long$variable <- NULL

# 合理重命名变量
colnames(stage_celltype_proportion_long) <- c("Developmental_Stage", "Cell_Type" , "Proportion")

# 绘制堆叠柱状图
stage_celltype_proportion_plot <- ggplot(
  stage_celltype_proportion_long,
  aes(x = Cell_Type, y = Proportion, fill = Developmental_Stage)
) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Proportion", fill = "Developmental Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate axis labels for better readability
  RotatedAxis()

# 保存图片
ggsave(
  "stage_celltype_proportion_plot.png",
  plot = stage_celltype_proportion_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Display the plot
print(stage_celltype_proportion_plot)


# 生成细胞类型/发育时间点的计数表
celltype_stage_count <- t(stage_celltype_count)

# 计算每个发育时间点的总细胞数
stage_totals <- colSums(celltype_stage_count)

# 计算每种细胞类型在每个发育时间点中的比例
celltype_stage_proportion <- as.data.frame(t(t(celltype_stage_count) / stage_totals))

# 将数据转换为长格式
celltype_stage_proportion_long <- melt(celltype_stage_proportion)

celltype_stage_proportion_long$variable <- NULL

# 合理重命名变量
colnames(celltype_stage_proportion_long) <- c("Cell_Type", "Developmental_Stage", "Proportion")

# 绘制堆叠柱状图
celltype_stage_proportion_plot <- ggplot(
  celltype_stage_proportion_long,
  aes(x = Developmental_Stage, y = Proportion, fill = Cell_Type)
) +
  geom_bar(stat = "identity") +
  labs(x = "Developmental Stage", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate axis labels for better readability
  RotatedAxis()

# 保存图片
ggsave(
  "celltype_stage_proportion_plot.png",
  plot = celltype_stage_proportion_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Display the plot
print(celltype_stage_proportion_plot)
