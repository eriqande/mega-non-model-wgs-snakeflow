rule samtools_stats:
    input:
        get_bams_for_samtools_stats,
    output:
        "results/bqsr-round-{bqsr_round}/qc/samtools_stats/{sample}.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/samtools_stats/{sample}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/samtools_stats/{sample}.bmk",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input} > {output} 2> {log} "


# this is a version that can create the same output for bqsr_round > 0
# from the bam in the recal directory.  We'll see if this works or not...
# rule samtools_stats2:
#     input:
#         "results/bqsr-round-{bqsr_round}/recal/{sample}.bam",
#     output:
#         "results/bqsr-round-{bqsr_round}/qc/samtools_stats/{sample}.txt",
#     log:
#         "results/bqsr-round-{bqsr_round}/logs/samtools_stats/{sample}.log",
#     benchmark:
#         "results/bqsr-round-{bqsr_round}/benchmarks/samtools_stats/{sample}.bmk",
#     wrapper:
#         "v1.1.0/bio/samtools/stats"


rule multiqc_dir:
    input:
        get_multiqc_inputs
    output:
        "results/bqsr-round-{bqsr_round}/qc/multiqc.html",
        directory("results/bqsr-round-{bqsr_round}/qc/multiqc_data")
    log:
        "results/bqsr-round-{bqsr_round}/logs/multiqc.log",
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/multiqc/multiqc.bmk",
    resources:
        mem_mb = 36800
    wrapper:
        "v3.5.1-2-g5f72db4/bio/multiqc"



# here is a rule to use bcftools stats to summarize what is found in the
# final BCF files.  Basically, we want to run bcftools stats with the option
# to compile statistics about every single sample, and also to make a histogram
# of our NMISS INFO tag that we made earlier.  The wildcarding here is so that
# we can do this three different ways:  1) with all the variants 2) with only those
# variants that PASS the filter, and 3) with only those variants that don't pass
# the hard filtering. Note that fitering in the bcftools stats command seems to
# conflict with samples in there so that things don't work, so I am piping it
# THIS is done on each chromo or scaff group since doing it on the full file can
# take upwards of 5 hours on a big data set.
rule bcf_section_summaries:
    input:
        "results/bqsr-round-{bqsr_round}/hard_filtering/both-filtered-{sg_or_chrom}.bcf",
    output:
        "results/bqsr-round-{bqsr_round}/qc/bcftools_stats/sections/{sg_or_chrom}-{filter_condition}.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/bcftools_stats/{sg_or_chrom}-{filter_condition}.log",
    params:
        comma_samples=",".join(unique_sample_ids),
        filter_opt=get_bcftools_stats_filter_option,
        stop=len(unique_sample_ids),
        steps=len(unique_sample_ids) + 1
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/bcftools_stats/{sg_or_chrom}-{filter_condition}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view {params.filter_opt} -Ou {input} | "
        " bcftools stats -s {params.comma_samples} "
        " -u NMISS:0:{params.stop}:{params.steps} > "
        " {output} 2> {log} "



# this is messed up.  Too much duplication and bullshit, but I
# want to do the maf-filtered ones, too.
rule bcf_maf_section_summaries:
    input:
        "results/bqsr-round-{bqsr_round}/hard_filtering/both-filtered-{sg_or_chrom}-maf-{maf}.bcf",
    output:
        "results/bqsr-round-{bqsr_round}/qc/bcftools_stats/maf_sections/{sg_or_chrom}-maf-{maf}.txt",
    log:
        "results/bqsr-round-{bqsr_round}/logs/bcftools_stats/{sg_or_chrom}-maf-{maf}.log",
    params:
        comma_samples=",".join(unique_sample_ids),
        stop=len(unique_sample_ids),
        steps=len(unique_sample_ids) + 1
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/bcftools_stats/{sg_or_chrom}-maf-{maf}.bmk",
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools stats -s {params.comma_samples} "
        " -u NMISS:0:{params.stop}:{params.steps}  {input} > "
        " {output} 2> {log} "


rule combine_bcftools_stats:
    input:
        expand("results/bqsr-round-{{bqsr_round}}/qc/bcftools_stats/sections/{sgc}-{{filter_condition}}.txt", sgc=unique_chromosomes + unique_scaff_groups)
    output:
        "results/bqsr-round-{bqsr_round}/qc/bcftools_stats/all-{filter_condition}.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/bqsr-round-{bqsr_round}/logs/combine_bcftools_stats/{filter_condition}.log",
    shell:
        " plot-vcfstats -m {input} > {output} 2>{log} "


rule combine_maf_bcftools_stats:
    input:
        expand("results/bqsr-round-{{bqsr_round}}/qc/bcftools_stats/maf_sections/{sgc}-maf-{{maf}}.txt", sgc=unique_chromosomes + unique_scaff_groups)
    output:
        "results/bqsr-round-{bqsr_round}/qc/bcftools_stats/all-pass-maf-{maf}.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/bqsr-round-{bqsr_round}/logs/combine_maf_bcftools_stats/maf-{maf}.log",
    shell:
        " plot-vcfstats -m {input} > {output} 2>{log} "
