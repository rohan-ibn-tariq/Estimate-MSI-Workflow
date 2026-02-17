# ============================================
# Microsatellite Detection (Dynamic WXS/WGS)
# ============================================


rule pytrf_find_repeats:
    """
Find all tandem repeats (microsatellites) in reference genome.

Detects mono- through hexanucleotide repeats.
Runs once for entire genome (used by both WXS and WGS).

Matches MSIsensor-pro scan defaults:
    -l 8: mononucleotide min 8 bases  -> mono=8
    -r 5: di-hexa min 5 repeats       -> di/tri/tetra/penta/hexa=5
    -s 6: max motif length = 6        -> findstr covers mono-hexa

Parameters:
    min_repeats: [8,5,5,5,5,5]
    - Mononucleotide: 8 repeats (e.g., AAAAAAAA)
    - Di/tri/tetra/penta/hexa: 5 repeats each
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
    threads: config["threads"]["single"]
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


# rule pytrf_filter_homopolymers:
#     """
# Filter homopolymers exceeding maximum length threshold.

# Matches MSIsensor-pro scan -m 50 (maximal homopolymer size).

# PyTRF CSV columns (NO header):
#     1. Chromosome
#     2. Start (1-based)
#     3. End (1-based)
#     4. Motif
#     5. Motif length
#     6. Repeat number
#     7. Repeat length

# Filter logic:
#     - For mononucleotides (motif_length == 1):
#         Keep only if repeat_length <= max_homopolymer
#     - For di-hexa (motif_length >= 2):
#         Keep all (no maximum defined in MSIsensor-pro)
# """
#     input:
#         csv="resources/microsatellites/all_repeats.csv",
#     output:
#         csv="resources/microsatellites/all_repeats_filtered.csv",
#     params:
#         max_homopolymer=config["msi"]["max_homopolymer"],
#     log:
#         "logs/microsatellites/pytrf_filter.log",
#     threads: config["threads"]["single"]
#     shell:
#         """
#         awk -F',' '
#             $5 == 1 && $7 > {params.max_homopolymer} {{
#                 next
#             }}
#             {{
#                 print
#             }}
#         ' {input.csv} \
#         > {output.csv} \
#         2> {log}
#         """


rule pytrf_to_bed:
    """
Convert PyTRF CSV output to UCSC BED Scheme
format (excluding bin column) and sort.

Note: PyTRF CSV has NO header line.
Columns in PyTRF CSV:
1. Chromosome
2. Start (1-based, inclusive)
3. End (1-based, inclusive)
4. Motif
5. Motif length
6. Repeat Number
7. Repeat Length
"""
    input:
        csv="resources/microsatellites/all_repeats.csv",
    output:
        bed=temp("resources/microsatellites/all_microsatellites.unsorted.bed"),
    log:
        "logs/microsatellites/pytrf_to_bed.log",
    threads: config["threads"]["single"]
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


rule extract_exonic_regions:
    """
Extract exonic coordinates from GTF for WXS filtering.

1. Filters GTF for lines where feature type (column 3) = "exon"
2. Extracts: chromosome (col 1), start (col 4 - 1), end (col 5)
3. Sorts lexicographically using bedtools

Example GTF line:
1  ensembl  exon  1000  2000  .  ...

Example BED output:
1  999  2000

Note: GTF is 1-based(start & end both inclusive),
      BED is 0-based start (start inclusive, end exclusive).
"""
    input:
        gtf="resources/genome/genome.gtf",
    output:
        bed="resources/microsatellites/exons.bed",
    log:
        "logs/microsatellites/extract_exons.log",
    conda:
        "../envs/bedtools.yaml"
    threads: config["threads"]["single"]
    shell:
        """
        awk '$3 == "exon" {{print $1 "\t" $4-1 "\t" $5}}' {input.gtf} | \
        bedtools sort -i - \
        > {output.bed} \
        2> {log}
        """


rule create_wxs_microsatellites:
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
        left="resources/microsatellites/all_microsatellites.bed",
        right="resources/microsatellites/exons.bed",
    output:
        "resources/microsatellites/wxs_microsatellites.bed",
    params:
        extra="-wa -u",
    log:
        "logs/microsatellites/create_wxs.log",
    wrapper:
        "v8.1.1/bio/bedtools/intersect"


rule create_wgs_microsatellites:
    """
Create WGS (genome-wide) microsatellite BED.

Simply copies all microsatellites (already sorted).
Separate rule for clarity and dynamic workflow.
"""
    input:
        bed="resources/microsatellites/all_microsatellites.bed",
    output:
        bed="resources/microsatellites/wgs_microsatellites.bed",
    log:
        "logs/microsatellites/create_wgs.log",
    threads: config["threads"]["single"]
    shell:
        """
cp {input.bed} {output.bed} 2> {log}
"""
