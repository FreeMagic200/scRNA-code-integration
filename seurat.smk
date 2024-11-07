# Seurat规则定义
rule seurat_read_data:
    input:
        rds = config["input"]["seurat"]
    output:
        rds = config["processed_dir"]["seurat"] + "/01_read_data.rds"
    log:
        "logs/seurat/01_read_data.log"
    params:
        input_file = config["input"]["seurat"],
        output_file = config["processed_dir"]["seurat"] + "/01_read_data.rds"
    script:
        "workflow/scripts/seurat/01_read_data.R"

rule seurat_normalize_stat:
    input:
        rds = config["processed_dir"]["seurat"] + "/01_read_data.rds"
    output:
        rds = config["processed_dir"]["seurat"] + "/02_normalized_stat.rds"
    params:
        figures_dir = config["figures_dir"]["seurat"]
    log:
        "logs/seurat/02_normalized_stat.log"
    script:
        "workflow/scripts/seurat/02_normalize_stat.R"

rule seurat_qc_filter:
    input:
        rds = config["processed_dir"]["seurat"] + "/02_normalized_stat.rds"
    output:
        rds = config["processed_dir"]["seurat"] + "/03_filtered_data.rds"
    params:
        qc = config["qc"],
    log:
        "logs/seurat/03_qc_filter.log"
    script:
        "workflow/scripts/seurat/03_qc_filter.R"

rule seurat_doublet_finder:
    input:
        rds = config["processed_dir"]["seurat"] + "/03_filtered_data.rds"
    output:
        rds = config["processed_dir"]["seurat"] + "/04_doublet_data.rds"
    params:
        figures_dir = config["figures_dir"]["seurat"],
        n_var_features = config["n_var_features"],
        n_pcs = config["n_pcs"],
        percent_doublets = config["qc"]["percent_doublets"],
        regress_features = config["regress_features"]
    log:
        "logs/seurat/04_doublet_finder.log"
    script:
        "workflow/scripts/seurat/04_doublet_finder.R"

rule seurat_postQC_process:
    input:
        rds = config["processed_dir"]["seurat"] + "/04_doublet_data.rds"
    output:
        rds = config["processed_dir"]["seurat"] + "/05_postQC_data.rds"
    params:
        n_var_features = config["n_var_features"],
        n_pcs = config["n_pcs"],
        regress_features = config["regress_features"]
    log:
        "logs/seurat/05_postQC_process.log"
    script:
        "workflow/scripts/seurat/05_postQC_process.R"

rule seurat_integrate:
    input:
        rds = config["processed_dir"]["seurat"] + "/05_postQC_data.rds"
    output:
        rds = config["processed_dir"]["seurat"] + "/06_integrated_data.rds"
    log:
        "logs/seurat/06_integrate.log"
    script:
        "workflow/scripts/seurat/06_integrate.R"
