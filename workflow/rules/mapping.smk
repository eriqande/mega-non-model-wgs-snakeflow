
# eca commented out the trimlog from this as it is huge and
# not particularly useful by default, as far as I can tell.
rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}---{unit}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}---{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{sample}---{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{sample}---{unit}.2.unpaired.fastq.gz"),
        #trimlog="results/trimmed/{sample}---{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        #extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    benchmark:
        "results/benchmarks/trim_reads_pe/{sample}---{unit}.bmk"
    log:
        "results/logs/trim_reads_pe/{sample}---{unit}.log",
    wrapper:
        "v1.1.0/bio/trimmomatic/pe"


# eca modified this.  The idea is to give 4 threads to bwa.
# and it will get 4 cores and also take all the memory you'd
# expect for those cores.  Sedna's machines are almost all
# 20 core units, so this should fill them up OK.
rule map_reads:
    input:
        reads = [
            "results/trimmed/{sample}---{unit}.1.fastq.gz",
            "results/trimmed/{sample}---{unit}.2.fastq.gz"
        ],
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}---{unit}.sorted.bam"),
    log:
        "results/logs/map_reads/{sample}---{unit}.log",
    benchmark:
        "results/benchmarks/map_reads/{sample}---{unit}.bmk"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    resources:
        time = "23:59:59",
    threads: 4
    wrapper:
        "0.59.2/bio/bwa/mem"



rule mark_duplicates:
    input:
        get_all_bams_of_common_sample
    output:
        bam=protected("results/mkdup/{sample}.bam"),
        bai=protected("results/mkdup/{sample}.bai"),
        metrics="results/qc/mkdup/{sample}.metrics.txt",
    log:
        "results/logs/picard/mkdup/{sample}.log",
    benchmark:
        "results/benchmarks/mark_duplicates/{sample}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1
    wrapper:
        "v1.1.0/bio/picard/markduplicates"


