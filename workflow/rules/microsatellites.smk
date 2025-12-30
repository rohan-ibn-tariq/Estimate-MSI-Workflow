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
        csv="resources/microsatellites/all_repeats.csv",
    params:
        min_repeats=config["microsatellites"]["min_repeats"],
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
