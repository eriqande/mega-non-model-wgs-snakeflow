
rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.fastq.gz"),
        r2=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.fastq.gz"),
        #r1_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.unpaired.fastq.gz"),
        #r2_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.unpaired.fastq.gz"),
        html="results/bqsr-round-{bqsr_round}/qc/fastp/{sample}---{unit}.html",
        json="results/bqsr-round-{bqsr_round}/qc/fastp/{sample}---{unit}.json"
    conda:
        "../envs/fastp.yaml"
    log:
        out="results/bqsr-round-{bqsr_round}/logs/trim_reads_pe/{sample}---{unit}.log",
        err="results/bqsr-round-{bqsr_round}/logs/trim_reads_pe/{sample}---{unit}.err"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/trim_reads_pe/{sample}---{unit}.bmk"
    params:
        trim_settings=config["params"]["fastp"]["pe"]["trimmer"],
    shell:
        " fastp -i {input.r1} -I {input.r2} "
        "       -o {output.r1} -O {output.r2} "
        "       -h {output.html} -j {output.json} "
        "  {params.trim_settings} > {log.out} 2> {log.err} "


# eca modified this.  The idea is to give 4 threads to bwa.
# and it will get 4 cores and also take all the memory you'd
# expect for those cores.  Sedna's machines are almost all
# 20 core units, so this should fill them up OK.
rule map_reads:
    input:
        reads = [
            "results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.fastq.gz",
            "results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.fastq.gz"
        ],
        idx=rules.bwa_index.output,
    output:
        temp("results/bqsr-round-{bqsr_round}/mapped/{sample}---{unit}.sorted.bam"),
    log:
        "results/bqsr-round-{bqsr_round}/logs/map_reads/{sample}---{unit}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/map_reads/{sample}---{unit}.bmk"
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
        sort_extra=""
    threads: 4
    resources:
        mem_mb=19200,
        time="23:59:59"
    wrapper:
        "v1.23.3/bio/bwa/mem"



rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam",
        bai="results/bqsr-round-{bqsr_round}/mkdup/{sample}.bai",
        metrics="results/bqsr-round-{bqsr_round}/qc/mkdup/{sample}.metrics.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/picard/mkdup/{sample}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1
    wrapper:
        "v1.1.0/bio/picard/markduplicates"


