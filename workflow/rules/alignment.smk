# ============================================
# Read Alignment and BAM Processing
# ============================================


def get_read_group(wildcards):
    """
    Generate read group string from sample metadata.

    Read groups track sample origin and sequencing details.
    Required for proper BAM file headers.

    Returns:
        str: BWA-formatted read group string
    """
    sample_row = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]

    rg_id = sample_row["ena_experiment_accession"]  # SRX11633435
    rg_sm = wildcards.sample  # Colo205
    rg_pl = sample_row["platform"]  # ILLUMINA
    rg_pm = sample_row["instrument_model"]  # Illumina NovaSeq 6000

    return rf"-R '@RG\tID:{rg_id}\tSM:{rg_sm}\tPL:{rg_pl}\tPM:{rg_pm}'"


rule bwa_mem:
    """
Align paired-end reads to reference genome using BWA-MEM2.
Sorting is done during alignment.
"""
    input:
        reads=[
            "resources/samples/{sample}_R1.fastq.gz",
            "resources/samples/{sample}_R2.fastq.gz",
        ],
        idx=multiext(
            "resources/genome/genome.fasta",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
            ".0123",
        ),
    output:
        "results/mapped/{sample}.bam",
    params:
        extra=lambda wildcards: get_read_group(wildcards),
        sort="samtools",
        sort_order="coordinate",
    log:
        "logs/alignment/{sample}_bwa.log",
    threads: config["threads"]["alignment"]
    wrapper:
        "v8.1.1/bio/bwa-mem2/mem"


rule samtools_index:
    """
Index BAM file for random access.
Enables fast access to specific genomic regions.
"""
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/mapped/{sample}.bam.bai",
    log:
        "logs/alignment/{sample}_index.log",
    threads: config["threads"]["default"]
    wrapper:
        "v8.1.1/bio/samtools/index"
