rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "results/logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v1.1.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "results/logs/genome-faidx.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools faidx {input}"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "results/logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "




rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/logs/bwa_index.log",
    resources:
        mem_mb=36900,
    cache: True
    wrapper:
        "0.59.2/bio/bwa/index"

