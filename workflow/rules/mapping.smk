rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        "results/logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.59.2/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="results/trimmed/{sample}-{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        "results/logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.59.2/bio/trimmomatic/pe"


# eca modified this.  The idea is to give 4 threads to bwa.
# and it will get 4 cores and also take all the memory you'd
# expect for those cores.  Sedna's machines are almost all
# 20 core units, so this should fill them up OK.
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam"),
    log:
        "results/logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    resources:
        cpus = 4,
        mem_mb = 18400  
    threads: 4
    wrapper:
        "0.59.2/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=protected("results/mkdup/{sample}-{unit}.bam"),
        metrics="results/qc/mkdup/{sample}-{unit}.metrics.txt",
    log:
        "results/logs/picard/mkdup/{sample}-{unit}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
        java_opts = "-Xmx4G ",
    resources:
        cpus = 1,
        mem_mb = 4600
    wrapper:
        "0.59.2/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        bam=protected("recal/{sample}-{unit}.bam"),
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    wrapper:
        "0.59.2/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "results/logs/samtools/index/{prefix}.log",
    wrapper:
        "0.59.2/bio/samtools/index"
