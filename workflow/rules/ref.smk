rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "results/bqsr-round-0/logs/get_genome.log",
    benchmark:
        "results/bqsr-round-0/benchmarks/get_genome/get_genome.bmk",
    params:
        url=config["ref"]["genome_url"],
    conda:
        "../envs/wget.yaml"
    shell:
        " (tmp_dir=$(mktemp -d) && "
        " URL={params.url} && "
        " if [[ $URL =~ \.gz$ ]]; then EXT='.gz'; else EXT=''; fi && "
        " wget -O $tmp_dir/file$EXT $URL && "
        " if [[ $URL =~ \.gz$ ]]; then gunzip $tmp_dir/file$EXT; fi && "
        " mv $tmp_dir/file {output}) > {log} 2>&1 "

rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "results/bqsr-round-0/logs/genome_faidx.log",
    benchmark:
        "results/bqsr-round-0/benchmarks/genome_faidx/genome_faidx.bmk",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "results/bqsr-round-0/logs/genome_dict.log",
    benchmark:
        "results/bqsr-round-0/benchmarks/genome_dict/genome_dict.bmk"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools dict {input} > {output} 2> {log} "




rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/bqsr-round-0/logs/bwa_index.log",
    benchmark:
        "results/bqsr-round-0/benchmarks/bwa_index/bwa_index.bmk",
    resources:
        mem_mb=36900,
    wrapper:
        "0.59.2/bio/bwa/index"


rule bowtie2_index:
    input:
        fasta="resources/genome.fasta",
    params:
        prefix="resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log:
        "results/bqsr-round-0/logs/bowtie2_index.log",
    benchmark:
        "results/bqsr-round-0/benchmarks/bowtie2_index/bowtie2_index.bmk",
    resources:
        mem_mb=36900,
    shell:
        """
        bowtie2-build {input.fasta} {params.prefix} > {log} 2>&1
        """


