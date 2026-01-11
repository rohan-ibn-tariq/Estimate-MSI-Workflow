# ============================================
# MS Coverage Analysis
# ============================================


rule sort_wxs_bed_to_match_genome:
    """
    Sort WXS microsatellite BED to match genome order.
    """
    input:
        bed="resources/microsatellites/wxs_microsatellites.bed",
        genome="resources/genome/genome.fasta.fai",
    output:
        bed="resources/coverage/wxs_microsatellites.sorted.bed",
    log:
        "logs/coverage/sort_wxs_bed.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools sort -faidx {input.genome} -i {input.bed} > {output.bed} 2> {log}
        """


rule sort_wgs_bed_to_match_genome:
    """
    Sort WGS microsatellite BED to match genome order.
    """
    input:
        bed="resources/microsatellites/wgs_microsatellites.bed",
        genome="resources/genome/genome.fasta.fai",
    output:
        bed="resources/coverage/wgs_microsatellites.sorted.bed",
    log:
        "logs/coverage/sort_wgs_bed.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools sort -faidx {input.genome} -i {input.bed} > {output.bed} 2> {log}
        """


rule compute_ms_coverage:
    """
Compute coverage fraction for each MS region.

Output columns (bedtools coverage default):
1-4: chr, start, end, name (from BED)
5: overlap_count (number of reads overlapping)
6: bases_covered (number of bases with non-zero coverage)
7: length (total length of MS region)
8: fraction (bases_covered / length)
"""
    input:
        bed=lambda wildcards: get_sorted_microsatellite_bed(wildcards),
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        genome="resources/genome/genome.fasta.fai",
    output:
        coverage="results/coverage/{sample}/ms_coverage.tsv",
    log:
        "logs/coverage/{sample}_coverage.log",
    conda:
        "../envs/bedtools.yaml"
    threads: config["threads"]["default"]
    resources:
        mem_mb=100
    shell:
        """
        bedtools coverage \
            -sorted \
            -g {input.genome} \
            -a {input.bed} \
            -b {input.bam} \
        > {output.coverage} \
        2> {log}
        """


rule plot_coverage_histogram:
    """
Plot histogram of MS region coverage fractions.
"""
    input:
        coverage="results/coverage/{sample}/ms_coverage.tsv",
    output:
        plot="results/coverage/{sample}/coverage_histogram.html",
    params:
        title=lambda wildcards: f"MS Coverage Distribution - {wildcards.sample}",
    log:
        "logs/coverage/{sample}_plot_coverage.log",
    conda:
        "../envs/coverage.yaml"
    script:
        "../scripts/plot_coverage_histogram.py"
