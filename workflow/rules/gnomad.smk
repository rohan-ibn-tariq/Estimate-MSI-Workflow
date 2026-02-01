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


rule annotate_sample_with_gnomad:
    """Annotate sample VCF with gnomAD population frequencies"""
    input:
        vcf="results/calls/{sample}.vcf.gz",
        csi="results/calls/{sample}.vcf.gz.csi"
    output:
        vcf="results/calls/{sample}_annotated.vcf.gz",
        csi="results/calls/{sample}_annotated.vcf.gz.csi"
    params:
        pop=lambda wildcards: get_gnomad_population(wildcards),
        gender=lambda wildcards: get_sample_gender(wildcards),
        gnomad_type=lambda wildcards: "exomes" if get_process_type(wildcards) == "WXS" else "genomes",
        s3_prefix="s3://prefix-path/resources/gnomad", # Replace with actual S3 path
        workdir="dev_sandbox/gnomad_annotation_{sample}"
    log:
        "logs/gnomad/annotate_{sample}.log"
    conda:
        "../envs/gnomad_annotation.yaml"
    shell:
        """
        # Determine sex-specific suffix for X/Y chromosomes
        if [ "{params.gender}" = "male" ]; then
            SEX_SUFFIX="_XY"
        elif [ "{params.gender}" = "female" ]; then
            SEX_SUFFIX="_XX"
        else
            SEX_SUFFIX=""
        fi

        echo "Starting gnomAD annotation for {wildcards.sample}..." > {log}
        echo "Population: {params.pop}" >> {log}
        echo "Gender: {params.gender}" >> {log}
        echo "Type: {params.gnomad_type}" >> {log}

        # Create clean working directory
        mkdir -p {params.workdir}
        cd {params.workdir}

        # Copy original to working file
        cp ../../{input.vcf} working.vcf.gz
        cp ../../{input.csi} working.vcf.gz.csi

        # Annotate autosomes (1-22) with general population AF
        for chr in {{1..22}}; do
            echo "Annotating chr$chr..." >> ../../{log}

            # Download gnomAD file for this chromosome
            s5cmd cp {params.s3_prefix}/{params.gnomad_type}/gnomad.{params.gnomad_type}.v4.1.sites.${{chr}}.renamed.vcf.gz gnomad_chr${{chr}}.vcf.gz 2>> ../../{log}
            s5cmd cp {params.s3_prefix}/{params.gnomad_type}/gnomad.{params.gnomad_type}.v4.1.sites.${{chr}}.renamed.vcf.gz.csi gnomad_chr${{chr}}.vcf.gz.csi 2>> ../../{log}

            # Annotate (only variants on this chromosome)
            bcftools annotate \
              -r $chr \
              -a gnomad_chr${{chr}}.vcf.gz \
              -c INFO/POPULATION_AF:=INFO/{params.pop} \
              -h <(echo '##INFO=<ID=POPULATION_AF,Number=A,Type=Float,Description="gnomAD population allele frequency">') \
              working.vcf.gz \
              -Oz -o working_new.vcf.gz 2>> ../../{log}

            # Replace working file
            rm working.vcf.gz working.vcf.gz.csi
            mv working_new.vcf.gz working.vcf.gz
            bcftools index working.vcf.gz 2>> ../../{log}

            # Delete gnomAD file
            rm gnomad_chr${{chr}}.vcf.gz gnomad_chr${{chr}}.vcf.gz.csi
        done

        # Annotate sex chromosomes + MT
        for chr_info in "MT:genomes/gnomad.genomes.v3.1.sites.MT.renamed.vcf.gz:{params.pop}" \
                        "X:{params.gnomad_type}/gnomad.{params.gnomad_type}.v4.1.sites.X.renamed.vcf.gz:{params.pop}${{SEX_SUFFIX}}" \
                        "Y:{params.gnomad_type}/gnomad.{params.gnomad_type}.v4.1.sites.Y.renamed.vcf.gz:{params.pop}${{SEX_SUFFIX}}"; do

            chr=$(echo $chr_info | cut -d: -f1)
            gnomad_path=$(echo $chr_info | cut -d: -f2)
            af_field=$(echo $chr_info | cut -d: -f3)

            echo "Annotating $chr..." >> ../../{log}

            # Download gnomAD file
            s5cmd cp {params.s3_prefix}/${{gnomad_path}} gnomad_${{chr}}.vcf.gz 2>> ../../{log}
            s5cmd cp {params.s3_prefix}/${{gnomad_path}}.csi gnomad_${{chr}}.vcf.gz.csi 2>> ../../{log}

            # Annotate
            bcftools annotate \
              -r $chr \
              -a gnomad_${{chr}}.vcf.gz \
              -c INFO/POPULATION_AF:=INFO/${{af_field}} \
              working.vcf.gz \
              -Oz -o working_new.vcf.gz 2>> ../../{log}

            # Replace working file
            rm working.vcf.gz working.vcf.gz.csi
            mv working_new.vcf.gz working.vcf.gz
            bcftools index working.vcf.gz 2>> ../../{log}

            # Delete gnomAD file
            rm gnomad_${{chr}}.vcf.gz gnomad_${{chr}}.vcf.gz.csi
        done

        # Move final output to destination
        mv working.vcf.gz ../../{output.vcf}
        mv working.vcf.gz.csi ../../{output.csi}

        # Clean up working directory
        cd ../..
        rm -rf {params.workdir}

        echo "Annotation complete!" >> {log}
        """
