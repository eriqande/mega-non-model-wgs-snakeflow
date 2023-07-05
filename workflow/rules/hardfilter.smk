

# this is the hard-filtering routine recommended by
# the folks that develop the GATK.

# It marks sites as PASS in the FILTER column,
# or it writes out the name of the particular filters
# that were not passed.

rule make_snp_vcf:
    input:
        vcf="results/bqsr-round-{bqsr_round}/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/logs/gatk/selectvariants/select-snps-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/make_snp_vcf/selectvariants-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "



rule make_indel_vcf:
    input:
        vcf="results/bqsr-round-{bqsr_round}/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/bqsr-round-{bqsr_round}/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/logs/gatk/selectvariants/select-indels-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/make_indel_vcf/selectvariants-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type INDEL -O {output.vcf} > {log} 2>&1 "




rule hard_filter_snps:
    input:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/snps-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/snps-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/logs/gatk/variantfiltration/snps-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/hard_filter_snps/variantfiltration-snps-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        "  -filter 'QD < 2.0' --filter-name 'QD2' "
        "  -filter 'QUAL < 30.0' --filter-name 'QUAL30' "
        "  -filter 'SOR > 3.0' --filter-name 'SOR3' "
        "  -filter 'FS > 60.0' --filter-name 'FS60' "
        "  -filter 'MQ < 40.0' --filter-name 'MQ40' "
        "  -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' "
        "  -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' "
        " -O {output.vcf} > {log} 2>&1 "




rule hard_filter_indels:
    input:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/indels-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/indels-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        idx="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/bqsr-round-{bqsr_round}/logs/gatk/variantfiltration/indels-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/hard_filter_indels/variantfiltration-indels-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        "gatk VariantFiltration "
        " -V {input.vcf} "
        "  -filter 'QD < 2.0' --filter-name 'QD2' "
        "  -filter 'QUAL < 30.0' --filter-name 'QUAL30' "
        "  -filter 'FS > 200.0' --filter-name 'FS200' "
        "  -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' "
        " -O {output.vcf} > {log} 2>&1 "




rule bung_filtered_vcfs_back_together:
    input:
        snp="results/bqsr-round-{bqsr_round}/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz",
        indel="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        snp_idx="results/bqsr-round-{bqsr_round}/hard_filtering/snps-filtered-{sg_or_chrom}.vcf.gz.tbi",
        indel_idx="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi"
    output:
        vcf="results/bqsr-round-{bqsr_round}/hard_filtering/both-filtered-{sg_or_chrom}.bcf",
    log:
        "results/bqsr-round-{bqsr_round}/logs/bung_filtered_vcfs_back_together/bung-{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/bung_filtered_vcfs_back_together/bcftools-{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools concat -a {input.snp} {input.indel} | "
        " bcftools view -Ob > {output.vcf}; ) 2> {log} "


rule maf_filter:
    input:
        "results/bqsr-round-{bqsr_round}/hard_filtering/both-filtered-{sg_or_chrom}.bcf"
    output:
        "results/bqsr-round-{bqsr_round}/hard_filtering/both-filtered-{sg_or_chrom}-maf-{maf}.bcf"
    log:
        "results/bqsr-round-{bqsr_round}/logs/maf_filter/{sg_or_chrom}-maf-{maf}.log",
    params:
        maf="{maf}"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/maf_filter/{sg_or_chrom}-maf-{maf}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools view -Ob -i 'FILTER=\"PASS\" & MAF > {params.maf} ' "
        " {input} > {output} 2>{log} "




