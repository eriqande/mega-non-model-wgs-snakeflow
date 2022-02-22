rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "results/logs/get_genome.log",
    benchmark:
        "results/benchmarks/get_genome/get_genome.bmk",
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
        "results/logs/genome_faidx.log",
    benchmark:
        "results/benchmarks/genome_faidx/genome_faidx.bmk",
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
        "results/logs/genome_dict.log",
    benchmark:
        "results/benchmarks/genome_dict/genome_dict.bmk"
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
    benchmark:
        "results/benchmarks/bwa_index/bwa_index.bmk",
    resources:
        mem_mb=36900,
    cache: True
    wrapper:
        "0.59.2/bio/bwa/index"

