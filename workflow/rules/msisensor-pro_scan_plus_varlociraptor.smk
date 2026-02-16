# ============================================
# MSIsensor-pro scan + Varlociraptor MSI detection
# ============================================


rule msisensor_pro_scan:
    """
Run MSIsensor-pro scan to identify microsatellite sites in the reference genome.
    """
    input:
        reference="resources/genome/genome.fasta",
    output:
        msisensor_sites="resources/msisensor-pro/microsatellites/msisensor_pro_sites.list",
    log:
        "logs/msisensor-pro/msisensor_pro_scan.log",
    conda:
        "../envs/msisensor-pro.yaml"
    shell:
        """
        msisensor-pro scan -d {input.reference} -o {output.msisensor_sites}  2>> {log}
        """


rule convert_msisensor_to_bed:
    """
    Convert MSIsensor-pro sites to 0-based BED format, calculate 
    the end position, format the name as 'repeatnumxMotif', 
    and sort using bedtools.
    """
    input:
        msisensor_sites="resources/msisensor-pro/microsatellites/msisensor_pro_sites.list"
    output:
        msisensor_bed="resources/msisensor-pro/microsatellites/wgs_msisensor_microsatellites.bed"
    log:
        "logs/msisensor-pro/convert_to_bed.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        # awk fields explanation:
        # $1 = chromosome
        # Note: MSIsensor-pro output is 0-based and also BED format is 0-based, so we keep $2 as is for chromStart.
        # $2 = chromStart (0-based)
        # $2+($3*$4) = chromEnd (exclusive)
        # $4"x"$5 = name (repeatnumxMotif)
        
        awk 'NR > 1 {{print $1 "\t" $2 "\t" ($2 + ($3 * $4)) "\t" $4 "x" $5}}' {input.msisensor_sites} \
        | bedtools sort -i - > {output.msisensor_bed} 2> {log}
        """


rule create_msisensor_wxs_microsatellites:
    """
Create WXS (exome) microsatellite BED.

1. Takes all microsatellites (genome-wide)
2. Keeps only those that overlap exonic regions

Note: Output is already sorted (bedtools intersect preserves order)

Flags:
-wa: Write original entry from A (keep full microsatellite info)
-u: Write each A entry only once (unique)
"""
    input:
        left="resources/msisensor-pro/microsatellites/wgs_msisensor_microsatellites.bed",
        right="resources/microsatellites/exons.bed",
    output:
        "resources/msisensor-pro/microsatellites/wxs_msisensor_microsatellites.bed",
    params:
        extra="-wa -u",
    log:
        "logs/msisensor-pro/create_wxs.log",
    wrapper:
        "v8.1.1/bio/bedtools/intersect"


rule estimate_msi_msisensor_pro_plus_varlociraptor:
    """
Estimate microsatellite instability using Varlociraptor.

Requirements:
- VCF: Lexicographically sorted
- BED: Lexicographically sorted
- Chromosomes: Matching naming (no "chr")
- Dynamic: Uses WXS or WGS BED based on sample type
"""
    input:
        vcf="results/calls/{sample}_annotated.vcf.gz",
        bed=lambda wildcards: get_microsatellite_msisensor(wildcards),
    output:
        pseudotime_data="results/msisensor-pro/msi/{sample}/pseudotime.tsv",
        pseudotime_plot="results/msisensor-pro/msi/{sample}/pseudotime.vl.json",
        distribution_data="results/msisensor-pro/msi/{sample}/distribution.tsv",
        distribution_plot="results/msisensor-pro/msi/{sample}/distribution.vl.json",
    log:
        "logs/msisensor-pro/{sample}_estimate.log",
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
