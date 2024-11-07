import os
import scanpy as sc
import pandas as pd
import numpy as np
from snakemake.logging import logger
from scipy.sparse import issparse

# 打印当前工作目录
current_working_directory = os.getcwd()
logger.info(f"当前工作目录: {current_working_directory}")

# 打印所有传递的命令行参数
input_file = snakemake.input[0]
output_file = snakemake.output[0]
qc_params = snakemake.config["qc"]
figures_dir = snakemake.config["figures_dir"]["scanpy"]  # 假设figures_dir是一个字典

logger.info(f"传递的参数:\n输入文件路径: {input_file}\n输出文件路径: {output_file}")
logger.info(f"QC 参数: {qc_params}")
logger.info(f"figures_dir: {figures_dir}")

try:
    # 从 H5AD 文件读取 AnnData 对象
    adata = sc.read_h5ad(input_file)
    logger.info(f"成功从以下路径读取 AnnData 对象: {input_file}")
    
    # 确认数据读取
    logger.info(f"AnnData 对象包含 {adata.n_obs} 行和 {adata.n_vars} 列的矩阵")
    logger.info(f"特征（变量）名称：{adata.var_names[:5]}")
    
    # 设置线粒体基因、血红蛋白基因和核糖体基因
    adata.var['mito'] = adata.var_names.str.startswith('MT-')
    adata.var['hb'] = adata.var_names.str.contains(r'^HB[^(egf)|^(s1l)|^(p1)].+')
    adata.var['ribo'] = adata.var_names.str.startswith('RP')
    
    # 计算 QC 指标
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mito', 'hb', 'ribo'], 
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # 确认计算结果
    logger.info(f"数据中的列名: {adata.obs.columns.tolist()}")
    
    # QC 过滤函数
    def qc_filter(adata, qc_params):
        logger.info("应用 QC 过滤条件")

        # 检查是否包含 Malat1 表达量信息
        if 'Malat1' not in adata.var_names:
            raise ValueError("Malat1 表达量信息缺失")

        # 分别计算每个条件
        cond1 = adata.obs['n_genes_by_counts'].values < qc_params['max_nFeature_RNA']
        cond2 = adata.obs['n_genes_by_counts'].values > qc_params['min_nFeature_RNA']
        cond3 = adata.obs['total_counts'].values < qc_params['max_nCount_RNA']
        cond4 = adata.obs['total_counts'].values > qc_params['min_nCount_RNA']
        cond5 = adata.obs['pct_counts_mito'].values < qc_params['max_percent_mt']
        cond6 = adata.obs['pct_counts_ribo'].values < qc_params['max_percent_ribo']
        cond7 = adata.obs['pct_counts_hb'].values < qc_params['max_percent_hb']
        
        # 处理 Malat1 表达量
        malat1_data = adata[:, 'Malat1'].X
        if issparse(malat1_data):
            malat1_data = malat1_data.toarray()
        cond8 = malat1_data.flatten() > qc_params['min_num_Malat1']
        
        # 将所有条件组合起来
        combined_conditions = cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8

        # 应用过滤条件
        adata_filtered = adata[combined_conditions]

        logger.info(f"过滤后的数据包含 {adata_filtered.n_obs} 行和 {adata_filtered.n_vars} 列")
        return adata_filtered
    
    # 应用 QC 过滤
    adata_filtered = qc_filter(adata, qc_params)
    
    # 保存过滤后的 AnnData 对象
    adata_filtered.write_h5ad(output_file)
    logger.info(f"过滤后的 AnnData 对象保存到: {output_file}")

except Exception as e:
    logger.error(f"发生错误: {str(e)}")
    raise

