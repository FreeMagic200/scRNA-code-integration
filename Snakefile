configfile: "config/config.yaml"
rule all:
    input:
        seurat_data = config["processed_dir"]["seurat"] + "/01_read_data.rds",
        scanpy_data = config["processed_dir"]["scanpy"] + "/01_read_data.h5ad"
    shell: 
      "echo Seurat data: {input.seurat_data}"
      "echo Scanpy data: {input.scanpy_data}"
# 包含 Seurat 和 Scanpy 的规则文件
include: "seurat.smk"
include: "scanpy.smk"
