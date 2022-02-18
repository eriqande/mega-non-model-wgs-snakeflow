rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}---{unit}.html",
        zip="results/qc/fastqc/{sample}---{unit}_fastqc.zip",
    log:
        "results/logs/fastqc/{sample}---{unit}.log",
    wrapper:
        "0.59.2/bio/fastqc"


rule samtools_stats:
    input:
        "results/bams_sampmerged/{sample}.bam",
    output:
        "results/qc/samtools-stats/{sample}.txt",
    log:
        "results/logs/samtools-stats/{sample}.log",
    wrapper:
        "0.59.2/bio/samtools/stats"


rule multiqc_and_index_bams:
    input:
        expand("results/qc/samtools-stats/{sample}.txt", sample=sample_list),
        expand("results/qc/fastqc/{u.sample}---{u.unit}_fastqc.zip", u=units.itertuples()),
        list(dict.fromkeys(expand("results/qc/mkdup/{u.sample}---{u.library}.metrics.txt", u=units.itertuples()))),
        expand("results/bams_sampmerged/{sample}.bam.bai", sample=sample_list)  # only here to force indexing of bams in this early step
    output:
        report(
            "results/qc/multiqc.html",
            caption="../../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "results/logs/multiqc.log",
    wrapper:
        "0.59.2/bio/multiqc"
