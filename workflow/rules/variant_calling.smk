# ============================================
# Variant Calling with Varlociraptor
# ============================================

# wildcard_constraints:
#     sample="[^/.]+"


rule bcftools_mpileup:
    """
Generate pileup from BAM alignments.
First step of candidate variant calling.
"""
    input:
        alignments="results/mapped/{sample}.bam",
        ref="resources/genome/genome.fasta",
        index="resources/genome/genome.fasta.fai",
    output:
        "results/pileups/{sample}.pileup.bcf",
    params:
        uncompressed_bcf=False,
        extra="--max-depth 200 --min-BQ 20 -a AD,DP",
    log:
        "logs/bcftools/{sample}_mpileup.log",
    threads: config["threads"]["default"]
    wrapper:
        "v8.1.1/bio/bcftools/mpileup"


rule bcftools_call:
    """
Call candidate variants from pileup.
"""
    input:
        pileup="results/pileups/{sample}.pileup.bcf",
    output:
        calls="results/candidates/{sample}.vcf.gz",
    params:
        uncompressed_bcf=False,
        caller="-m",
        extra="--ploidy 2 --variants-only",
    log:
        "logs/bcftools/{sample}_call.log",
    threads: config["threads"]["default"]
    wrapper:
        "v8.1.1/bio/bcftools/call"


rule varlociraptor_alignment_properties:
    """
Estimate alignment properties for Varlociraptor.
"""
    input:
        ref="resources/genome/genome.fasta",
        ref_idx="resources/genome/genome.fasta.fai",
        alignments="results/mapped/{sample}.bam",
        aln_idx="results/mapped/{sample}.bam.bai",
    output:
        "results/alignment-properties/{sample}.json",
    log:
        "logs/varlociraptor/{sample}_alignment_properties.log",
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/varlociraptor/estimate-alignment-properties"


rule varlociraptor_preprocess:
    """
Preprocess alignments at candidate variant sites.
"""
    input:
        ref="resources/genome/genome.fasta",
        alignment_properties="results/alignment-properties/{sample}.json",
        alignments="results/mapped/{sample}.bam",
        candidate_variants="results/candidates/{sample}.vcf.gz",
    output:
        "results/observations/{sample}.bcf",
    log:
        "logs/varlociraptor/{sample}_preprocess.log",
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/varlociraptor/preprocess-variants"


rule varlociraptor_call:
    """
Call variants using Varlociraptor Bayesian model.
"""
    input:
        observations="results/observations/{sample}.bcf",
        scenario=config["varlociraptor"]["scenario"],
    output:
        "results/calls/{sample}.bcf",
    log:
        "logs/varlociraptor/{sample}_call.log",
    params:
        samples=["sample"],
        extra="",
    wrapper:
        "v8.1.1/bio/varlociraptor/call-variants"


rule bcftools_view_to_vcf:
    """
Convert BCF to VCF (required for bedtools).
"""
    input:
        "results/calls/{sample}.bcf",
    output:
        "results/calls/{sample}_tmp.vcf",
    log:
        "logs/bcftools/{sample}_view_vcf.log",
    params:
        extra="-O v",
    wrapper:
        "v8.1.1/bio/bcftools/view"


rule bedtools_sort_vcf:
    """
Sort VCF lexicographically and preserve header.
"""
    input:
        in_file="results/calls/{sample}_tmp.vcf",
    output:
        "results/calls/{sample}.vcf",
    log:
        "logs/bedtools/{sample}_sort_vcf.log",
    params:
        extra="-header",
    wrapper:
        "v8.1.1/bio/bedtools/sort"
