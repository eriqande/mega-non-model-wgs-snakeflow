rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}_fastqc.zip",
    log:
        "results/logs/fastqc/{sample}-{unit}.log",
    wrapper:
        "0.59.2/bio/fastqc"


rule samtools_stats:
    input:
        "results/mkdup/{sample}-{unit}.bam",
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt",
    log:
        "results/logs/samtools-stats/{sample}-{unit}.log",
    wrapper:
        "0.59.2/bio/samtools/stats"


rule multiqc_and_index_bams:
    input:
        expand("results/qc/samtools-stats/{u.sample}-{u.unit}.txt", u=units.itertuples()),
        expand("results/qc/fastqc/{u.sample}-{u.unit}_fastqc.zip", u=units.itertuples()),
        expand("results/qc/mkdup/{u.sample}-{u.unit}.metrics.txt", u=units.itertuples()),
        expand("results/mkdup/{u.sample}-{u.unit}.bam.bai", u=units.itertuples())
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
