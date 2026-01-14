# ============================================
# MSIsensor-pro rules (MSI detection using MSIsensor-pro)
# ============================================


# def get_msisensor_sites(wildcards):
#     process_type = get_process_type(wildcards)
#     if process_type == "WXS":
#         return "resources/msisensor/microsatellites/msisensor_wxs.list"
#     else:
#         return "resources/msisensor/microsatellites/msisensor_wgs.list"


rule download_msisensor2_models:
    """
Download pretrained MSIsensor2 tumor-only models (hg38).
"""
    output:
        model_dir=directory("resources/msisensor2/models_hg38"),
    log:
        "logs/msisensor2/download_models.log",
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.model_dir})
        mkdir -p $(dirname {log})

        tmp_log=$(mktemp)
        tmp_dir=$(mktemp -d)

        git clone --depth 1 https://github.com/niu-lab/msisensor2.git $tmp_dir >> $tmp_log 2>&1

        cd $tmp_dir
        git sparse-checkout init --cone >> $tmp_log 2>&1
        git sparse-checkout set models_hg38 >> $tmp_log 2>&1
        cd -

        mkdir -p {output.model_dir}
        cp -r "$tmp_dir"/models_hg38/* {output.model_dir}/ >> "$tmp_log" 2>&1

        ls -lh "{output.model_dir}" >> "$tmp_log" 2>&1
        echo "Download complete" >> "$tmp_log"

        rm -rf $tmp_dir

        cat "$tmp_log" > "{log}"
        rm "$tmp_log"
        """


rule msisensor2_scan:
    """
Scan reference for microsatellites matching PyTRF [5,3,3,3,3,3].

-l 5: Minimal homopolymer (mono) = 5
-s 6: Maximal MS unit = 6 (to include hexa like PyTRF)
-r 3: Minimal repeats = 3
# NOTE:
# This (full command) approximates PyTRF [5,3,3,3,3,3] as closely as possible
"""
    input:
        ref="resources/genome/genome.fasta",
    output:
        sites="resources/msisensor2/microsatellites_msisensor2.list",
    log:
        "logs/msisensor2/scan.log",
    conda:
        "../envs/msisensor.yaml"
    shell:
        """
        msisensor2 scan \
            -d {input.ref} \
            -o {output.sites} \
            -l 5 \
            -s 6 \
            -r 3 \
        > {log} 2>&1
        """


rule extract_exonic_regions_msisensor2:
    """
Extract exonic coordinates from GTF.
Output (4-column BED): chr  start  end  gene_name
"""
    input:
        gtf="resources/genome/genome.gtf",
    output:
        bed="resources/msisensor2/exons.bed",
    log:
        "logs/msisensor2/extract_exons.log",
    conda:
        "../envs/bedtools.yaml"
    threads: config["threads"]["single"]
    shell:
        """
        awk '$3 == "exon" {{
            match($9, /gene_name "([^"]+)"/, arr);
            gene = (arr[1] != "" ? arr[1] : "unknown");
            print $1 "\t" ($4-1) "\t" $5 "\t" gene
        }}' {input.gtf} | \
        bedtools sort -i - \
        > {output.bed} \
        2> {log}
        """


rule msisensor2_msi:
    """
MSIsensor2 tumor-only MSI detection.

For WXS: Uses exon BED to restrict analysis
For WGS: Analyzes whole genome (no -e flag)

Parameters matching PyTRF [5,3,3,3,3,3]:
-l 5: Minimal homopolymer
-p 5: Minimal homopolymer for distribution
-q 3: Minimal MS unit size
-s 3: Minimal MS size for distribution
-c 20: Coverage threshold (WXS: 20, WGS: 15)
"""
    input:
        bam="results/mapped/{sample}.bam",
        models="resources/msisensor2/models_hg38",
        exons="resources/msisensor2/exons.bed",  # Always required for dependency
    output:
        msi="results/msisensor2/{sample}",
        dis="results/msisensor2/{sample}_dis",
        somatic="results/msisensor2/{sample}_somatic",
    params:
        prefix=lambda wildcards: f"results/msisensor2/{wildcards.sample}",
        process_type=get_process_type,
        exon_flag=lambda wildcards: (
            "-e resources/msisensor2/exons.bed"
            if get_process_type(wildcards) == "WXS"
            else ""
        ),
        coverage_threshold=lambda wildcards: (
            20 if get_process_type(wildcards) == "WXS" else 15
        ),
    log:
        "logs/msisensor2/{sample}_msi.log",
    conda:
        "../envs/msisensor.yaml"
    threads: 12  # MSIsensor2 can multithread, @TODO: make configurable
    shell:
        """
        mkdir -p results/msisensor2

        msisensor2 msi \
            -M {input.models} \
            -t {input.bam} \
            {params.exon_flag} \
            -o {params.prefix} \
            -b {threads} \
            -c {params.coverage_threshold} \
            -l 5 \
            -p 5 \
            -q 3 \
            -s 3 \
        > {log} 2>&1
        """
