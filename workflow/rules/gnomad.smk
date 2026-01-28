# ============================================
# Gnomad (Dynamic WXS/WGS)
# ============================================


# GNOMAD_VERSION_MAIN = "4.1"
# GNOMAD_VERSION_MT = "3.1"
CHROMOSOMES_MAIN = [str(i) for i in range(1, 23)] + ["X", "Y"]

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
        tbi="resources/gnomad/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz.tbi",
    params:
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{gnomad_type}/gnomad.{gnomad_type}.v4.1.sites.chr{chr}.vcf.bgz",
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        wget -q -O {output.vcf} {params.url}
        wget -q -O {output.tbi} {params.url}.tbi
        """


rule download_gnomad_mt:
    output:
        vcf="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
        tbi="resources/gnomad/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi",
    params:
        url="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        wget -q -O {output.vcf} {params.url}
        wget -q -O {output.tbi} {params.url}.tbi
        """
