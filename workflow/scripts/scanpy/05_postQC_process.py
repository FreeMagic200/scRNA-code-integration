import os
import scanpy as sc
import matplotlib.pyplot as plt
from snakemake.logging import logger

# 打印当前工作目录
current_working_directory = os.getcwd()
logger.info(f"当前工作目录: {current_working_directory}")

# 获取输入和输出文件路径
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# 获取参数
n_var_features = snakemake.config["n_var_features"]
n_pcs = snakemake.config["n_pcs"]
figures_dir = snakemake.config["figures_dir"]["scanpy"]  # 获取 scanpy 图像保存目录
regress_features = snakemake.config["regress_features"]

logger.info(f"输入文件: {input_file}")
logger.info(f"输出文件: {output_file}")
logger.info(f"选择的高变基因数量: {n_var_features}")
logger.info(f"PCA 成分数量: {n_pcs}")
logger.info(f"回归特征: {regress_features}")

# 打印 figures_dir 以检查路径
logger.info(f"图像保存目录: {figures_dir}")

# 特征名称转换表
feature_name_mapping = {
    "nCount_RNA": "total_counts",
    "nFeature_RNA": "n_genes_by_counts",
    "percent.mt": "pct_counts_mito",
    "percent.ribo": "pct_counts_ribo",
    "percent.hb": "pct_counts_hb",
    "G2M_score": "G2M_score",
    "S_score": "S_score"
    # 如果有其他特征需要转换，请在这里添加
}

# 转换特征名称
regress_features = [feature_name_mapping.get(f, f) for f in regress_features]
logger.info(f"转换后的回归特征: {regress_features}")

# 读取细胞周期基因文件
with open('mmus_s_genes.txt') as f:
    s_genes = [line.strip() for line in f]
with open('mmus_g2m_genes.txt') as f:
    g2m_genes = [line.strip() for line in f]

logger.info(f"读取的 S 期基因: {s_genes}")
logger.info(f"读取的 G2M 期基因: {g2m_genes}")

# 确保 figures_dir 目录存在
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)
    logger.info(f"创建目录: {figures_dir}")

try:
    # 读取 AnnData 对象
    adata = sc.read_h5ad(input_file)
    logger.info(f"成功读取 AnnData 对象，包含 {adata.n_obs} 个细胞和 {adata.n_vars} 个基因")

    # 过滤掉 doublets
    adata = adata[~adata.obs['predicted_doublet'], :]
    logger.info(f"过滤 doublets 后，剩余 {adata.n_obs} 个细胞")

    # 细胞周期打分
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    logger.info("细胞周期打分已完成")

    # 寻找高变基因
    sc.pp.highly_variable_genes(adata, n_top_genes=n_var_features)
    logger.info(f"已选择 {n_var_features} 个高变基因")

    # 确保图像保存目录存在
#    figures_path = os.path.join(figures_dir)
 #   if not os.path.exists(figures_path):
  #      os.makedirs(figures_path)
   #     logger.info(f"创建目录: {figures_path}")

    # 绘制高变基因图
    #high_var_genes_file = os.path.join(figures_path, 'highly_variable_genes.png')
    #sc.pl.highly_variable_genes(adata, save=high_var_genes_file)
    #logger.info(f"高变基因图已保存至 {high_var_genes_file}")

    # 只保留高变基因
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    logger.info(f"只保留高变基因后，数据包含 {adata.n_vars} 个基因")

    # 数据标准化并进行回归
    sc.pp.regress_out(adata, regress_features)
    logger.info("数据已进行回归处理")
    sc.pp.scale(adata, max_value=10)
    logger.info("数据已标准化")

    # 运行 PCA
    sc.tl.pca(adata, n_comps=n_pcs)
    logger.info(f"已运行 PCA，保留 {n_pcs} 个主成分")

    # 绘制 PCA 方差比例图
    #pca_variance_ratio_file = os.path.join(figures_dir, 'pca_variance_ratio.png')
    #sc.pl.pca_variance_ratio(adata, log=True, save=pca_variance_ratio_file)
    #logger.info(f"PCA 方差比例图已保存至 {pca_variance_ratio_file}")

    # 保存结果
    adata.write_h5ad(output_file)
    logger.info(f"处理后的 AnnData 对象已保存至 {output_file}")

except Exception as e:
    logger.error(f"发生错误: {str(e)}")
    raise

