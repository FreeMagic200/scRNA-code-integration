import numpy as np
import scanpy as sc
import scvi
from typing import Tuple
from snakemake.logging import logger

scvi.settings.seed = 42

# 获取输入和输出文件路径
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# 获取参数
layer_num = snakemake.params.n_layers
latent_num = snakemake.params.n_latent
epochs_num = snakemake.params.n_epochs
use_cuda = snakemake.params.use_cuda
batch_index = snakemake.params.batch_idx

logger.info(f"输入文件: {input_file}")
logger.info(f"输出文件: {output_file}")
logger.info(f"scVI 潜在空间维度: {latent_num}")
logger.info(f"训练轮次: {epochs_num}")
logger.info(f"是否使用 GPU: {use_cuda}")

# 读取 AnnData 对象
adata = sc.read_h5ad(input_file)
logger.info(f"成功读取 AnnData 对象，包含 {adata.n_obs} 个细胞和 {adata.n_vars} 个基因")

# SCVI setup
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_index)

model = scvi.model.SCVI(adata, n_layers = layer_num, n_latent = latent_num, gene_likelihood = "nb")

model.train(accelerator = "gpu")

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

SCVI_MDE_KEY = "X_scVI_MDE"
adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])

# 保存结果
adata.write_h5ad(output_file)
logger.info(f"处理后的 AnnData 对象已保存至 {output_file}")

