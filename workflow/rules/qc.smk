rule fastqc_read1:
    input:
        lambda wc: (get_fastq(wc))["r1"],
    output:
        html="results/qc/fastqc/{sample}---{unit}_R1.html",
        zip="results/qc/fastqc/{sample}---{unit}_R1_fastqc.zip",
    log:
        "results/logs/fastqc/{sample}---{unit}_R1.log",
    benchmark:
        "results/benchmarks/fastqc_read1/{sample}---{unit}_R1.bmk",
    wrapper:
        "v1.1.0/bio/fastqc"


rule fastqc_read2:
    input:
        lambda wc: (get_fastq(wc))["r2"],
    output:
        html="results/qc/fastqc/{sample}---{unit}_R2.html",
        zip="results/qc/fastqc/{sample}---{unit}_R2_fastqc.zip",
    log:
        "results/logs/fastqc/{sample}---{unit}_R2.log",
    benchmark:
        "results/benchmarks/fastqc_read2/{sample}---{unit}_R2.bmk",
    wrapper:
        "v1.1.0/bio/fastqc"



rule samtools_stats:
    input:
        "results/mkdup/{sample}.bam",
    output:
        "results/qc/samtools_stats/{sample}.txt",
    log:
        "results/logs/samtools_stats/{sample}.log",
    benchmark:
        "results/benchmarks/samtools_stats/{sample}.bmk",
    wrapper:
        "v1.1.0/bio/samtools/stats"


rule multiqc:
    input:
        expand("results/qc/samtools_stats/{sample}.txt", sample=sample_list),
        expand("results/qc/fastqc/{u.sample}---{u.unit}_R{r}_fastqc.zip", u=units.itertuples(), r = [1, 2]),
        list(dict.fromkeys(expand("results/qc/mkdup/{u.sample}.metrics.txt", u=units.itertuples()))),
        expand("results/logs/trim_reads_pe/{u.sample}---{u.unit}.log", u=units.itertuples()),
    output:
        "results/qc/multiqc.html",
    log:
        "results/logs/multiqc.log",
    benchmark:
        "results/benchmarks/multiqc/multiqc.bmk",
    resources:
        mem_mb = 36800
    wrapper:
        "v1.3.1/bio/multiqc"



# here is a rule to use bcftools stats to summarize what is found in the
# final VCF files.  Basically, we want to run bcftools stats with the option
# to compile statistics about every single sample, and also to make a histogram
# of our NMISS INFO tag that we made earlier.  The wildcarding here is so that
# we can do this three different ways:  1) with all the variants 2) with only those
# variants that PASS the filter, and 3) with only those variants that don't pass
# the hard filtering. Note that fitering in the bcftools stats command seems to
# conflict with samples in there so that things don't work, so I am piping it
rule vcf_summaries:
    input:
        "results/vcf/all-filtered.vcf.gz",
    output:
        "results/qc/bcftools_stats/bcftools-stats-{filter_condition}.txt",
    log:
        "results/logs/bcftools_stats/{filter_condition}.log",
    params:
        comma_samples=",".join(unique_sample_ids),
        filter_opt=get_bcftools_stats_filter_option,
        stop=len(unique_sample_ids),
        steps=len(unique_sample_ids) + 1
    benchmark:
        "results/benchmarks/bcftools_stats/{filter_condition}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view {params.filter_opt} -Ou {input} | "
        " bcftools stats -s {params.comma_samples} "
        " -u NMISS:0:{params.stop}:{params.steps} > "
        " {output} 2> {log} "


