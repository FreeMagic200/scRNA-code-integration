rule scanpy_read_data:
    input:
        h5ad= config["input"]["scanpy"]
    output:
        h5ad= config["processed_dir"]["scanpy"] + "/01_read_data.h5ad"
    log:
        "logs/scanpy/01_read_data.log"
    params:
        input_file=config["input"]["scanpy"]
    script:
        "workflow/scripts/scanpy/01_read_data.py"

rule scanpy_normalize_data:
    input:
        h5ad = config["processed_dir"]["scanpy"] + "/01_read_data.h5ad"
    output:
        h5ad = config["processed_dir"]["scanpy"] + "/02_normalized_data.h5ad"
    log:
        "logs/scanpy/02_normalize_data.log"
    params:
        input_file=config["processed_dir"]["scanpy"] + "/01_read_data.h5ad"
    script:
        "workflow/scripts/scanpy/02_normalize_data.py"

rule scanpy_qc_filter:
    input:
        h5ad = config["processed_dir"]["scanpy"] + "/02_normalized_data.h5ad"
    output:
        h5ad = config["processed_dir"]["scanpy"] + "/03_qc_filtered_data.h5ad"
    log:
        "logs/scanpy/03_qc_filter.log"
    params:
        qc=config["qc"],
        figures_dir=config["figures_dir"]["scanpy"]
    script:
        "workflow/scripts/scanpy/03_qc_filter.py"

rule scanpy_doublet_detection:
    input:
        h5ad = config["processed_dir"]["scanpy"] + "/03_qc_filtered_data.h5ad"
    output:
        h5ad = config["processed_dir"]["scanpy"] + "/04_doublet_detected_data.h5ad"
    log:
        "logs/scanpy/04_doublet_detection.log"
    params:
        batch_key = "Batch",
        figures_dir = config["figures_dir"]["scanpy"]
    script:
        "workflow/scripts/scanpy/04_doublet_detection.py"

rule scanpy_postQC_process:
    input:
        h5ad = config["processed_dir"]["scanpy"] + "/04_doublet_detected_data.h5ad"
    output:
        h5ad = config["processed_dir"]["scanpy"] + "/05_postQC_processed_data.h5ad"
    log:
        "logs/scanpy/05_postQC_process.log"
    params:
        figures_dir = config["figures_dir"]["scanpy"]
    script:
        "workflow/scripts/scanpy/05_postQC_process.py"

# Define the rule to integrate scVI
rule scanpy_integrate:
    input:
        input_file=config["processed_dir"]["scanpy"] + "/05_postQC_processed_data.h5ad"
  # 输入的 AnnData 文件路径
    output:
        output_file="data/06_integrated.h5ad"  # 保存结果的 AnnData 文件路径
    params:
        n_layers=2,
        n_latent=2,
        n_epochs=400,
        use_cuda=True,
        batch_idx="Batch"
    script:
        "workflow/scripts/scanpy/06_integrate.py"


