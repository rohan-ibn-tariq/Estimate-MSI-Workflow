# ============================================
# MS Coverage Analysis
# ============================================


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
        bed=lambda wildcards: get_microsatellite_bed(wildcards),
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
    output:
        coverage="results/coverage/{sample}/ms_coverage.tsv",
    log:
        "logs/coverage/{sample}_coverage.log",
    conda:
        "../envs/bedtools.yaml"
    threads: config["threads"]["default"]
    shell:
        """
        bedtools coverage \
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
