# ============================================
# Reference Genome and Indices
# ============================================


rule get_genome:
    """Download reference genome from Ensembl."""
    output:
        "resources/genome/genome.fasta",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        chromosome=config["ref"].get("chromosome", ""),
    log:
        "logs/ref/get_genome.log",
    wrapper:
        "v8.1.1/bio/reference/ensembl-sequence"


rule get_annotation:
    """Download genome annotation from Ensembl."""
    output:
        "resources/genome/genome.gtf",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    log:
        "logs/ref/get_annotation.log",
    wrapper:
        "v8.1.1/bio/reference/ensembl-annotation"


rule genome_faidx:
    """Create FASTA index."""
    input:
        "resources/genome/genome.fasta",
    output:
        "resources/genome/genome.fasta.fai",
    log:
        "logs/ref/genome_faidx.log",
    wrapper:
        "v8.1.1/bio/samtools/faidx"


rule bwa_index:
    """Create BWA index for alignment."""
    input:
        "resources/genome/genome.fasta",
    output:
        multiext("resources/genome/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/ref/bwa_index.log",
    threads: config["threads"]["default"]
    resources:
        mem_mb=15000,
    wrapper:
        "v8.1.1/bio/bwa/index"
