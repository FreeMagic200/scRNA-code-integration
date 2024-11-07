library(Seurat)
library(ggplot2)
library(ggridges)
library(scales)

# 获取当前脚本所在的路径
script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 设置工作目录
setwd(script_path)

#配色来源 Liu Y, Zhang Q, Xing B, Luo N, Gao R, Yu K, Hu X, Bu Z, Peng J, Ren X, Zhang Z. Immune phenotypic linkage between colorectal cancer and liver metastasis. Cancer Cell. 2022 Apr 11;40(4):424-437.e5. doi: 10.1016/j.ccell.2022.02.013. Epub 2022 Mar 17. PMID: 35303421.
my_cols <- c(
  "#dc8e97",
  "#e3d1db",
  "#74a893",
  "#ac9141",
  "#5ac6e9",
  "#ebce8e",
  "#e5c06e",
  "#7587b1",
  "#c7deef",
  "#e97371",
  "#e1a4c6",
  "#916ba6",
  "#cb8f82",
  "#7db3af",
  "#d2e0ac"
)

obj <- readRDS("../../../data/processed/seurat/02_normalized_stat.rds")

# 提取数据并转换为数据框
gene_data <- as.data.frame(FetchData(obj, vars = c("nFeature_RNA", "Batch"), layer = "data"))

# 绘制 ridge plot

gene_ridge_plot <- ggplot(gene_data, aes(x = nFeature_RNA, y = Batch, fill = Batch)) +
  stat_density_ridges(geom = "density_ridges_gradient", scale = 2.3) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_vline(xintercept = 500,
             linetype = "dashed",
             color = "blue") +
  geom_vline(xintercept = 8000,
             linetype = "dashed",
             color = "red") +
  labs(title = "Distribution of Gene Number", x = "Number of Genes", y = "Batch") +
  scale_x_log10(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black")
  ) +
  annotate(
    "text",
    x = 500,
    y = Inf,
    label = paste0("500"),
    vjust = 1,
    hjust = -0.1,
    color = "blue"
  ) +
  annotate(
    "text",
    x = 8000,
    y = Inf,
    label = paste0("8,000"),
    vjust = 1,
    hjust = -0.1,
    color = "red"
  )

print(gene_ridge_plot)

# 提取数据并转换为数据框
umi_data <- as.data.frame(FetchData(obj, vars = c("nCount_RNA", "Batch"), layer = "data"))

# 绘制 ridge plot

umi_ridge_plot <- ggplot(umi_data, aes(x = nCount_RNA, y = Batch, fill = Batch)) +
  stat_density_ridges(geom = "density_ridges_gradient", scale = 2.3) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_vline(xintercept = 1000,
             linetype = "dashed",
             color = "blue") +
  geom_vline(xintercept = 30000,
             linetype = "dashed",
             color = "red") +
  labs(title = "Distribution of Normalized UMI Counts", x = "UMI Count Number", y = "Batch") +
  scale_x_log10(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black") 
  ) +
  annotate(
    "text",
    x = 1000,
    y = Inf,
    label = paste0("1,000"),
    vjust = 1,
    hjust = -0.1,
    color = "blue"
  ) +
  annotate(
    "text",
    x = 30000,
    y = Inf,
    label = paste0("30,000"),
    vjust = 1,
    hjust = -0.1,
    color = "red"
  )


print(umi_ridge_plot)


# 提取数据并转换为数据框
malat1_data <- as.data.frame(FetchData(obj, vars = c("Malat1", "Batch"), layer = "data"))

# 绘制 ridge plot
malat1_ridge_plot <- ggplot(malat1_data, aes(x = Malat1, y = Batch, fill = Batch)) +
  stat_density_ridges(geom = "density_ridges_gradient", scale = 2.3) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_vline(xintercept = 3,
             linetype = "dashed",
             color = "blue") +
  labs(title = "Distribution of Malat1 Counts", x = "Malat1 Count Number", y = "Batch") +
  scale_x_sqrt(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black")
  ) +
  annotate(
    "text",
    x = 3,
    y = Inf,
    label = paste0("3"),
    vjust = 1,
    hjust = -1.5,
    color = "blue"
  )


print(malat1_ridge_plot)


# featch data
mt_data<-as.data.frame(FetchData(obj, vars = c("percent.mt", "Batch"), layer = "data"))

# Create violin plot
mt_box_plot <- ggplot(mt_data, aes(x = Batch, y = percent.mt, fill = Batch)) +
  geom_boxplot(alpha = 0.95, outlier.size = 0.25) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_hline(yintercept = 10,
             linetype = "dashed",
             color = "red") +
  labs(title = "Percentage of Mitochondrial Counts", x = "Batch", y = "Percentage of Mitochondrial Counts (%)") +
  scale_y_sqrt(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  annotate(
    "text",
    x = Inf,
    y = 10,
    label = paste0("10"),
    vjust = -0.5,
    hjust = 1.1,
    color = "red"
  )

print(mt_box_plot)

# featch data
ribo_data<-as.data.frame(FetchData(obj, vars = c("percent.ribo", "Batch"), layer = "data"))

# Create violin plot
ribo_box_plot <- ggplot(ribo_data, aes(x = Batch, y = percent.ribo, fill = Batch)) +
  geom_boxplot(alpha = 0.95, outlier.size = 0.25) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_hline(yintercept = 20,
             linetype = "dashed",
             color = "red") +
  labs(title = "Percentage of Ribosomal Counts", x = "Batch", y = "Percentage of Ribosomal Counts (%)") +
  scale_y_sqrt(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  annotate(
    "text",
    x = Inf,
    y = 20,
    label = paste0("20"),
    vjust = -0.5,
    hjust = 1.1,
    color = "red"
  )

print(ribo_box_plot)

# featch data
hb_data<-as.data.frame(FetchData(obj, vars = c("percent.hb", "Batch"), layer = "data"))

# Create violin plot
hb_box_plot <- ggplot(hb_data, aes(x = Batch, y = percent.hb, fill = Batch)) +
  geom_boxplot(alpha = 0.95, outlier.size = 0.25) +
  scale_fill_manual(values = alpha(my_cols, 0.95), name = "Batch") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "red") +
  labs(title = "Percentage of Hemoglobin Counts", x = "Batch", y = "Percentage of Hemoglobin Counts (%)") +
  # scale_y_sqrt(labels = scales::comma) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(face = "bold", color = "grey"),
    plot.caption = element_text(color = "grey"),
    panel.background = element_blank(),
    # Remove panel background
    panel.grid.major = element_blank(),
    # Remove major grid lines
    panel.grid.minor = element_blank(),
    # Remove minor grid lines
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    # Retain axis ticks with black color
    axis.text = element_text(color = "black"),  # Retain axis text with black color
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1,
    label = paste0("1"),
    vjust = -0.5,
    hjust = 1,
    color = "red"
  )

print(hb_box_plot)

library(patchwork)

# 组合图表，两行三列
combined_plot <- (gene_ridge_plot + umi_ridge_plot + malat1_ridge_plot) /
  (mt_box_plot + ribo_box_plot + hb_box_plot) +
  plot_layout(nrow = 2)

# 保存组合后的图表为PDF
ggsave("combined_quality_control_plots.pdf", combined_plot, width = 24, height = 16, units = "in", device = "pdf")

obj$rpca_cell_type
