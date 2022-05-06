rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "results/logs/get_genome.log",
    benchmark:
        "results/benchmarks/get_genome/get_genome.bmk",
    params:
        url=config["ref"]["genome_url"],
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
        

rule define_chromosomes_and_scaffolds:
    input:
        fai="resources/genome.fasta.fai"
    output:
        chrom="resources/chromosomes.tsv",
        scaff="resources/scaffold_groups.tsv"
    log:
        "results/logs/define_chromosomes_and_scaffolds/define_chromosomes_and_scaffolds.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/make-chromosomes-and-scaffolds.R"

