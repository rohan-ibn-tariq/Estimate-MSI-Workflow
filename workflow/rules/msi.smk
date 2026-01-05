# ============================================
# MSI Estimation
# ============================================


rule estimate_msi:
    """
Estimate microsatellite instability using Varlociraptor.

Requirements:
- VCF: Lexicographically sorted
- BED: Lexicographically sorted
- Chromosomes: Matching naming (no "chr")
- Dynamic: Uses WXS or WGS BED based on sample type
"""
    input:
        vcf="results/calls/{sample}.vcf",
        bed=lambda wildcards: get_microsatellite_bed(wildcards),
    output:
        pseudotime_data="results/msi/{sample}/pseudotime.tsv",
        pseudotime_plot="results/msi/{sample}/pseudotime.vl.json",
        distribution_data="results/msi/{sample}/distribution.tsv",
        distribution_plot="results/msi/{sample}/distribution.vl.json",
    log:
        "logs/msi/{sample}_estimate.log",
    params:
        varlociraptor_bin=config["varlociraptor"]["binary"],
        varlociraptor_dir=config["varlociraptor"]["directory"],
        setup_script=config["varlociraptor"]["setup_script"],
        wrapper_script="workflow/scripts/varlociraptor-wrapper.sh",
        msi_threshold=config["msi"]["msi_threshold"],
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
