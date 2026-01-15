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
        awk -F'\t' '$3 == "exon" {{
            match($9, /gene_name[[:space:]]+"([^"]+)"/, arr);
            gene = (arr[1] != "" ? arr[1] : "unknown");
            print $1 "\t" ($4-1) "\t" $5 "\t" gene
        }}' {input.gtf} | \
        bedtools sort -i - \
        > {output.bed} \
        2> {log}
        """


rule download_chromToUcsc:
    """Download UCSC chromToUcsc Python script."""
    output:
        script="resources/msisensor2/chromToUcsc.py",
    log:
        "logs/msisensor2/download_chromToUcsc.log",
    shell:
        """
        wget -O {output.script} \
            https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc 2> {log}
        chmod +x {output.script}
        """


rule download_chromAlias:
    """Download hg38 chromAlias mapping file."""
    output:
        alias="resources/msisensor2/hg38.chromAlias.txt",
    log:
        "logs/msisensor2/download_chromAlias.log",
    shell:
        """
        wget -O {output.alias} \
            https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromAlias.txt 2> {log}
        """


rule convert_bam_to_ucsc:
    """Convert BAM chromosome names from Ensembl to UCSC format."""
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        script="resources/msisensor2/chromToUcsc.py",
        alias="resources/msisensor2/hg38.chromAlias.txt",
    output:
        bam="results/msisensor2/{sample}/converted.bam",
        bai="results/msisensor2/{sample}/converted.bam.bai",
    log:
        "logs/msisensor2/{sample}/convert_bam.log",
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]["indexing"]
    shell:
        """
        samtools view -h {input.bam} 2>> {log} | \
        python {input.script} -a {input.alias} -s 2>> {log} | \
        samtools view -b -@ {threads} -o {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """


rule convert_exon_bed_to_ucsc:
    """Convert exon BED chromosome names to UCSC format."""
    input:
        bed="resources/msisensor2/exons.bed",
        script="resources/msisensor2/chromToUcsc.py",
        alias="resources/msisensor2/hg38.chromAlias.txt",
    output:
        bed="resources/msisensor2/exon_ucsc.bed",
    log:
        "logs/msisensor2/convert_exon_bed.log",
    shell:
        """
        python {input.script} -a {input.alias} -i {input.bed} -o {output.bed} -s 2> {log}
        """


rule convert_msisensor2_list:
    """Convert MSIsensor2 list chromosome names to UCSC format."""
    input:
        list="resources/msisensor2/microsatellites_msisensor2.list",
        script="resources/msisensor2/chromToUcsc.py",
        alias="resources/msisensor2/hg38.chromAlias.txt",
    output:
        list="resources/msisensor2/microsatellites_list_ucsc.list",
    log:
        "logs/msisensor2/convert_list.log",
    shell:
        """
        python {input.script} -a {input.alias} -i {input.list} -o {output.list} -s 2> {log}
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
        bam="results/msisensor2/{sample}/converted.bam",
        bai="results/msisensor2/{sample}/converted.bam.bai",  # Needed for indexing
        models="resources/msisensor2/models_hg38",
        exons="resources/msisensor2/exon_ucsc.bed",  # Always required for dependency
        sites="resources/msisensor2/microsatellites_list_ucsc.list",
    output:
        msi="results/msisensor2/{sample}",
        dis="results/msisensor2/{sample}_dis",
        somatic="results/msisensor2/{sample}_somatic",
    params:
        prefix=lambda wildcards, output: output.msi,
        exon_flag=lambda wildcards, input: (
            f"-e {input.exons}" if get_process_type(wildcards) == "WXS" else ""
        ),
        coverage_threshold=lambda wildcards: (
            1 if get_process_type(wildcards) == "WXS" else 1
        ),
    log:
        "logs/msisensor2/{sample}_msi.log",
    conda:
        "../envs/msisensor.yaml"
    threads: 12  # MSIsensor2 can multithread, @TODO: make configurable
    shell:
        """
        mkdir -p $(dirname {output.msi})

        msisensor2 msi \
            -M {input.models} \
            -d {input.sites} \
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
