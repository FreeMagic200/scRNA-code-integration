import os
import scanpy as sc
import anndata as ad
from snakemake.logging import logger

# 打印当前工作目录
current_working_directory = os.getcwd()
logger.info(f"当前工作目录: {current_working_directory}")

# 打印所有传递的命令行参数
input_file = snakemake.input[0]
output_file = snakemake.output[0]
logger.info(f"传递的参数:\n输入文件路径: {input_file}\n输出文件路径: {output_file}")

try:
    # 从 H5AD 文件读取 AnnData 对象
    adata = sc.read_h5ad(input_file)
    logger.info(f"成功从以下路径读取 AnnData 对象: {input_file}")
    
    # 检查对象完整性
    if not isinstance(adata, ad.AnnData):
        raise ValueError("对象不是 'AnnData' 类")
    
    # 其他完整性检查示例
    if adata.X is None or adata.X.shape[0] == 0:
        raise ValueError("AnnData 对象中没有数据")
    if adata.obs is None or adata.obs.shape[0] == 0:
        raise ValueError("AnnData 对象中没有 metadata")
    
    logger.info("AnnData 对象通过完整性检查")
    
    # 保存处理后的 AnnData 对象
    adata.write_h5ad(output_file)
    logger.info(f"处理后的 AnnData 对象保存到: {output_file}")

except Exception as e:
    # 捕获错误并记录错误信息
    logger.error(f"发生错误: {str(e)}")
    raise

