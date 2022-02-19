
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
        "0.59.2/bio/trimmomatic/pe"


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
        cpus = 4,
        mem_mb = 18400  
    threads: 4
    wrapper:
        "0.59.2/bio/bwa/mem"


rule bam_merge_libraries:
    input: 
        bam=get_units_of_common_sample_and_lib
    output:
        bam=temp("results/bams_libmerged/{sample}---{library}.bam"),
        insize="results/benchmarks/bam_merge_libraries/{sample}---{library}.input_sizes"

    log:
        "results/logs/bam_merge_libraries/{sample}-{library}.log",
    benchmark:
        "results/benchmarks/bam_merge_libraries/{sample}---{library}.bmk"
    conda:
        "../envs/samtools.yaml"
    shell:
        "du -k {input} > {output.insize}; "
        "samtools merge {output.bam} {input.bam} 2> {log}"



rule bam_merge_samples:
    input: 
        bam=get_libmerged_bams_of_common_sample
    output:
        bam=protected("results/bams_sampmerged/{sample}.bam"),
        insize="results/benchmarks/bam_merge_samples/{sample}.input_sizes"
    log:
        "results/logs/bam_merge_samples/{sample}.log",
    benchmark:
        "results/benchmarks/bam_merge_samples/{sample}.bmk"
    conda:
        "../envs/samtools.yaml"
    shell:
        "du -k {input} > {output.insize}; "
        "samtools merge {output.bam} {input.bam} 2> {log}"



rule mark_duplicates:
    input:
        "results/bams_libmerged/{sample}---{library}.bam",
    output:
        bam=temp("results/mkdup/{sample}---{library}.bam"),
        metrics="results/qc/mkdup/{sample}---{library}.metrics.txt",
    log:
        "results/logs/picard/mkdup/{sample}---{library}.log",
    benchmark:
        "results/benchmarks/mark_duplicates/{sample}---{library}.bmk"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
    resources:
        cpus = 1,
        mem_mb = 4700
    wrapper:
        "v1.1.0/bio/picard/markduplicateswithmatecigar"







rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        bam=protected("recal/{sample}---{unit}.bam"),
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    log:
        "logs/gatk/bqsr/{sample}---{unit}.log",
    wrapper:
        "0.59.2/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "results/bams_sampmerged/{sample}.bam",
    output:
        protected("results/bams_sampmerged/{sample}.bam.bai"),
        insize="results/benchmarks/samtools_index/{sample}.input_sizes"
    log:
        "results/logs/samtools/index/{sample}.log",
    benchmark:
        "results/benchmarks/samtools_index/{sample}.bmk"
    conda:
        "../envs/samtools.yaml"
    shell:
        "du -k {input} > {output.insize}; "
        "samtools index {input} 2> {log}"
