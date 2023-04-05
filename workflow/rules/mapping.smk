
# eca commented out the trimlog from this as it is huge and
# not particularly useful by default, as far as I can tell.
rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.fastq.gz"),
        r2=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.2.unpaired.fastq.gz"),
        #trimlog="results/bqsr-round-{bqsr_round}/trimmed/{sample}---{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        #extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/trim_reads_pe/{sample}---{unit}.bmk"
    log:
        "results/bqsr-round-{bqsr_round}/logs/trim_reads_pe/{sample}---{unit}.log",
    wrapper:
        "v1.1.0/bio/trimmomatic/pe"


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
        "results/bqsr-round-{bqsr_round}/mapped/{sample}---{unit}.sorted.bam",
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
        time="10:00:00",
        mem_mb = 15200,
        cpus = 4,
	tmpdir = "results/bqsr-round-{bqsr_round}/mapped"
    wrapper:
        "v1.23.3/bio/bwa/mem"

rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam=protected("results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam"),
        bai=protected("results/bqsr-round-{bqsr_round}/mkdup/{sample}.bai"),
        metrics="results/bqsr-round-{bqsr_round}/qc/mkdup/{sample}.metrics.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/picard/mkdup/{sample}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        time="10:00:00",
	mem_mb = 38000,
        cpus = 10,
        tmpdir = "results/bqsr-round-0/tmp"
    wrapper:
        "v1.1.0/bio/picard/markduplicates"


