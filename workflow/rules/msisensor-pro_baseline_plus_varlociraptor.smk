# ============================================
# MSIsensor-pro baseline + Varlociraptor
# ============================================


rule download_msisensor_pro_baseline:
    """
Download official MSIsensor-pro baseline.
"""
    output:
        baseline="resources/msisensor-pro/baseline/GRCh38.baseline_TCGA-v1.3.tsv",
    log:
        "logs/msisensor-pro/download_baseline.log",
    shell:
        """
        wget -O {output.baseline} \
            https://raw.githubusercontent.com/xjtu-omics/msisensor-pro/refs/heads/master/data/GRCh38.baseline_TCGA-v1.3.tsv \
            2> {log}
        """


rule convert_baseline_to_bed:
    """
Convert MSIsensor-pro baseline TSV to sorted BED format.

Baseline TSV columns:
$1  = chromosome (chr1, chr2, ...)
$2  = location (1-based position)
$3  = repeat_unit_length
$4  = repeat_times
$5  = repeat_unit_bases

BED output (0-based coordinates):
chrom  chromStart  chromEnd  name
"""
    input:
        baseline="resources/msisensor-pro/baseline/GRCh38.baseline_TCGA-v1.3.tsv",
    output:
        bed="resources/msisensor-pro/baseline/wgs_microsatellites.bed",
    log:
        "logs/msisensor-pro/baseline_to_bed.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        awk -F'\t' 'NR > 1 {{
            gsub(/^chr/, "", $1)
            chrom = $1
            start = $2 - 1
            end = $2 + ($3 * $4) - 1
            name = $4 "x" $5
            print chrom "\t" start "\t" end "\t" name
        }}' {input.baseline} | \
        bedtools sort -i - > {output.bed} 2> {log}
        """


rule create_wxs_baseline_bed:
    """
Create WXS-specific baseline BED (only exonic sites).
"""
    input:
        left="resources/msisensor-pro/baseline/wgs_microsatellites.bed",
        right="resources/microsatellites/exons.bed",
    output:
        "resources/msisensor-pro/baseline/wxs_microsatellites.bed",
    params:
        extra="-wa -u",
    log:
        "logs/msisensor-pro/create_wxs_baseline.log",
    wrapper:
        "v8.1.1/bio/bedtools/intersect"


rule estimate_msi_with_msisensor_pro_baseline:
    """
Estimate MSI using Varlociraptor with MSIsensor-pro baseline.
"""
    input:
        vcf="results/calls/{sample}.vcf",
        bed=lambda wildcards: get_baseline_bed(wildcards),
    output:
        pseudotime_data="results/msisensor-pro-baseline-varlociraptor/msi/{sample}/pseudotime.tsv",
        pseudotime_plot="results/msisensor-pro-baseline-varlociraptor/msi/{sample}/pseudotime.vl.json",
        distribution_data="results/msisensor-pro-baseline-varlociraptor/msi/{sample}/distribution.tsv",
        distribution_plot="results/msisensor-pro-baseline-varlociraptor/msi/{sample}/distribution.vl.json",
    log:
        "logs/msisensor-pro-baseline-varlociraptor/{sample}_estimate.log",
    params:
        varlociraptor_bin=config["varlociraptor"]["binary"],
        varlociraptor_dir=config["varlociraptor"]["directory"],
        setup_script=config["varlociraptor"]["setup_script"],
        wrapper_script="workflow/scripts/varlociraptor-wrapper.sh",
        msi_threshold=config["msi"]["threshold"],
    threads: config["threads"]["varlociraptor"]
    shell:
        """
        mkdir -p $(dirname {output.pseudotime_data})
        
        bash {params.wrapper_script} \
            {params.varlociraptor_bin} \
            {params.varlociraptor_dir} \
            {params.setup_script} \
            estimate msi \
            {input.bed} \
            {input.vcf} \
            --threads {threads} \
            --msi-threshold {params.msi_threshold} \
            --data-pseudotime {output.pseudotime_data} \
            --plot-pseudotime {output.pseudotime_plot} \
            --data-distribution {output.distribution_data} \
            --plot-distribution {output.distribution_plot} \
            2> {log}
        """
