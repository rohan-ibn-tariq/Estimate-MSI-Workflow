# ============================================
# Variant Calling with Varlociraptor
# ============================================

wildcard_constraints:
    sample="[^/.]+"


rule freebayes_candidates:
    input:
        alns="results/mapped/{sample}.bam",
        idxs="results/mapped/{sample}.bam.bai",
        ref="resources/genome/genome.fasta",
        index="resources/genome/genome.fasta.fai",
    output:
        bcf=temp("results/candidates/{sample}.bcf")
    log:
        "logs/freebayes/{sample}.log"
    params:
        extra="--pooled-continuous -I -X"
    threads: config["threads"]["freebayes"]
    resources:
        mem_mb=16384
    wrapper:
        "v8.1.1/bio/freebayes"


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
        temp("results/alignment-properties/{sample}.json"),
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
        candidate_variants="results/candidates/{sample}.bcf",
    output:
        temp("results/observations/{sample}.bcf"),
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
        temp("results/calls/{sample}.bcf"),
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
        temp("results/calls/{sample}.tmp.vcf"),
    log:
        "logs/bcftools/{sample}_view_vcf.log",
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/bcftools/view"


rule bedtools_sort_vcf:
    """
Sort VCF lexicographically and preserve header.
"""
    input:
        in_file="results/calls/{sample}.tmp.vcf",
    output:
        "results/calls/{sample}.vcf",
    log:
        "logs/bedtools/{sample}_sort_vcf.log",
    params:
        extra="-header",
    wrapper:
        "v8.1.1/bio/bedtools/sort"
