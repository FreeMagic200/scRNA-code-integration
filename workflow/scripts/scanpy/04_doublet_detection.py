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
batch_key = snakemake.params.get("batch_key", "sample")
figures_dir = snakemake.config["figures_dir"]["scanpy"]

logger.info(f"输入文件: {input_file}")
logger.info(f"输出文件: {output_file}")
logger.info(f"Batch key: {batch_key}")

try:
    # 读取 AnnData 对象
    adata = sc.read_h5ad(input_file)
    logger.info(f"成功读取 AnnData 对象，包含 {adata.n_obs} 个细胞和 {adata.n_vars} 个基因")

    # 运行 Scrublet
    sc.pp.scrublet(adata, batch_key=batch_key)
    logger.info("完成 Scrublet 运行")

    # 打印 doublet 检测结果
    n_doublets = adata.obs['predicted_doublet'].sum()
    logger.info(f"检测到 {n_doublets} 个可能的 doublets")

    # 绘制 doublet score 分布图
    os.makedirs(figures_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(adata.obs['doublet_score'], bins=50)
    ax.set_xlabel('Doublet Score')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Doublet Scores')
    fig.savefig(os.path.join(figures_dir, 'doublet_score_distribution.png'))
    plt.close(fig)
    logger.info(f"Doublet score 分布图已保存至 {os.path.join(figures_dir, 'doublet_score_distribution.png')}")

    # 保存结果
    adata.write_h5ad(output_file)
    logger.info(f"已将包含 doublet 检测结果的 AnnData 对象保存至 {output_file}")

except Exception as e:
    logger.error(f"发生错误: {str(e)}")
    raise
