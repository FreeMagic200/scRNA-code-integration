snakemake --dag scanpy_integrate | dot -Tsvg > scanpy_dag.svg
snakemake --dag seurat_integrate | dot -Tsvg > seurat_dag.svg
