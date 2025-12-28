# ============================================
# Sample FASTQ Files Handling
# ============================================


def get_fastq_url(wildcards, read):
    """Get download URL for FASTQ from samples.tsv."""
    sample_row = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    return sample_row[f"url_{read}"]


rule download_fastq:
    """
Downloads paired-end FASTQ files from ENA using aria2c.
"""
    output:
        r1="resources/samples/{sample}_R1.fastq.gz",
        r2="resources/samples/{sample}_R2.fastq.gz",
    params:
        url_r1=lambda wildcards: get_fastq_url(wildcards, "r1"),
        url_r2=lambda wildcards: get_fastq_url(wildcards, "r2"),
        max_conn=4,
        split=4,
        min_split_size="1M",
        max_tries=5,
        retry_wait=3,
        timeout=60,
        connect_timeout=30,
        continue_download=True,
        allow_overwrite=True,
        auto_file_renaming=False,
        console_log_level="error",
        summary_interval=180,
    log:
        "logs/reads/{sample}_download.log",
    conda:
        "../envs/aria2.yaml"
    threads: 4
    retries: 3
    shell:
        """
        echo "Downloading {wildcards.sample} from ENA..." > {log}
        echo "R1: {params.url_r1}" >> {log}
        echo "R2: {params.url_r2}" >> {log}
        
        aria2c \
            --max-connection-per-server={params.max_conn} \
            --split={params.split} \
            --min-split-size={params.min_split_size} \
            --max-tries={params.max_tries} \
            --retry-wait={params.retry_wait} \
            --timeout={params.timeout} \
            --connect-timeout={params.connect_timeout} \
            --continue={params.continue_download} \
            --allow-overwrite={params.allow_overwrite} \
            --auto-file-renaming={params.auto_file_renaming} \
            --console-log-level={params.console_log_level} \
            --summary-interval={params.summary_interval} \
            --dir=$(dirname {output.r1}) \
            --out=$(basename {output.r1}) \
            {params.url_r1} \
            2>> {log} &
        
        aria2c \
            --max-connection-per-server={params.max_conn} \
            --split={params.split} \
            --min-split-size={params.min_split_size} \
            --max-tries={params.max_tries} \
            --retry-wait={params.retry_wait} \
            --timeout={params.timeout} \
            --connect-timeout={params.connect_timeout} \
            --continue={params.continue_download} \
            --allow-overwrite={params.allow_overwrite} \
            --auto-file-renaming={params.auto_file_renaming} \
            --console-log-level={params.console_log_level} \
            --summary-interval={params.summary_interval} \
            --dir=$(dirname {output.r2}) \
            --out=$(basename {output.r2}) \
            {params.url_r2} \
            2>> {log} &
        
        wait
        
        if [ ! -f {output.r1} ] || [ ! -f {output.r2} ]; then
            echo "ERROR: Download failed" >> {log}
            exit 1
        fi
        
        echo "Download complete" >> {log}
        echo "R1: $(du -h {output.r1} | cut -f1)" >> {log}
        echo "R2: $(du -h {output.r2} | cut -f1)" >> {log}
        """
