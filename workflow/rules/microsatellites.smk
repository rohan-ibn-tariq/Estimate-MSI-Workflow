# ============================================
# Microsatellite Detection (Dynamic WXS/WGS)
# ============================================


rule pytrf_find_repeats:
    """
Find all tandem repeats (microsatellites) in reference genome.

Detects mono- through hexanucleotide repeats.
Runs once for entire genome (used by both WXS and WGS).

Parameters:
    min_repeats: [5,3,3,3,3,3]
    - Mononucleotide: 5 repeats (e.g., AAAAA)
    - Di/tri/tetra/penta/hexa: 3 repeats each
"""
    input:
        fasta="resources/genome/genome.fasta",
    output:
        csv=protected("resources/microsatellites/all_repeats.csv"),
    params:
        min_repeats=config["msi"]["min_repeats"],
        fmt="csv",
    log:
        "logs/microsatellites/pytrf_find.log",
    conda:
        "../envs/pytrf.yaml"
    threads: config["threads"]["pytrf"]
    shell:
        """
        pytrf findstr \
            -r {params.min_repeats[0]} {params.min_repeats[1]} {params.min_repeats[2]} \
               {params.min_repeats[3]} {params.min_repeats[4]} {params.min_repeats[5]} \
            -f {params.fmt} \
            -o {output.csv} \
            {input.fasta} \
            > {log} 2>&1
        """


rule pytrf_to_bed:
    """
Convert PyTRF CSV output to UCSC BED Scheme
format (excluding bin column) and sort.

Note: PyTRF CSV has NO header line.
"""
    input:
        csv="resources/microsatellites/all_repeats.csv",
    output:
        bed=temp("resources/microsatellites/all_microsatellites.unsorted.bed"),
    log:
        "logs/microsatellites/pytrf_to_bed.log",
    threads: config["threads"]["pytrf"]
    shell:
        """
        awk -F',' '{{
            chrom = $1
            start = $2 - 1
            end = $3
            motif = $4
            copies = int($6)
            name = copies "x" motif
            print chrom "\t" start "\t" end "\t" name
        }}' {input.csv} \
        > {output.bed} \
        2> {log}
        """


rule bedtools_sort_all_microsatellites:
    """
Sort all microsatellites BED lexicographically.
Lexicographical order: 1, 10, 11, 2, 22, 3, X, Y
"""
    input:
        in_file="resources/microsatellites/all_microsatellites.unsorted.bed",
    output:
        "resources/microsatellites/all_microsatellites.bed",
    log:
        "logs/microsatellites/bedtools_sort_all.log",
    params:
        extra="",
    wrapper:
        "v8.1.1/bio/bedtools/sort"
