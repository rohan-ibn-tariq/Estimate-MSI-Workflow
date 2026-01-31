# ============================================
# Gnomad (Dynamic WXS/WGS)
# ============================================


# GNOMAD_VERSION_MAIN = "4.1"
# GNOMAD_VERSION_MT = "3.1"
CHROMOSOMES_MAIN = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHROMOSOMES_ALL = [str(i) for i in range(1, 23)] + ["MT", "X", "Y"]

def get_gnomad_types_needed():
    """Determine if we need exomes, genomes, or both"""
    active_samples = get_active_samples()
    process_types = samples_df[samples_df["sample"].isin(active_samples)][
        "process_as"
    ].unique()

    types = []
    if "WXS" in process_types:
        types.append("exomes")
    if "WGS" in process_types:
        types.append("genomes")

    return types


GNOMAD_TYPES = get_gnomad_types_needed()


rule download_gnomad_chr:
    output:
        vcf="resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz",
    params:
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz",
    log:
        "logs/gnomad/download_{gnomad_type}_chr{chr}.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        wget -q -O {output.vcf} {params.url} > {log} 2>&1
        """


rule download_gnomad_mt:
    output:
        vcf="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
    params:
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
    log:
        "logs/gnomad/download_genomes_chrM.log"
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        wget -q -O {output.vcf} {params.url} > {log} 2>&1
        """


rule download_all_gnomad:
    """Download all needed gnomAD files"""
    input:
        expand(
            "resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz",
            gnomad_type=GNOMAD_TYPES,
            chr=CHROMOSOMES_MAIN
        ),
        "resources/gnomad/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"


rule rename_gnomad_chr_main:
    input:
        vcf="resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz",
        rename_map="config/gnomad_chr_rename.txt"
    output:
        vcf="resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.{chr}.renamed.vcf.gz",
        csi="resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.{chr}.renamed.vcf.gz.csi"
    log:
        "logs/gnomad/rename_{gnomad_type}_chr{chr}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate --rename-chrs {input.rename_map} {input.vcf} -Oz -o {output.vcf} > {log} 2>&1
        bcftools index {output.vcf} >> {log} 2>&1
        """


rule rename_gnomad_mt:
    input:
        vcf="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
        rename_map="config/gnomad_chr_rename.txt"
    output:
        vcf="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.MT.renamed.vcf.gz",
        csi="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.MT.renamed.vcf.gz.csi"
    log:
        "logs/gnomad/rename_genomes_chrM.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate --rename-chrs {input.rename_map} {input.vcf} -Oz -o {output.vcf} > {log} 2>&1
        bcftools index {output.vcf} >> {log} 2>&1
        """
