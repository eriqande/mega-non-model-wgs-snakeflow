

# this is the hard-filtering routine recommended by
# the folks that develop the GATK.

# It marks sites as PASS in the FILTER column,
# or it writes out the name of the particular filters
# that were not passed.

rule make_snp_vcf:
    input:
        vcf="results/vcf/all.vcf.gz"
    output:
        vcf="results/hard_filtering/snps.vcf.gz"
    log:
        "results/logs/gatk/selectvariants/select-snps.log",
    benchmark:
        "results/benchmarks/make_snp_vcf/selectvariants-snps.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type SNP -O {output.vcf} > {log} 2>&1 "



rule make_indel_vcf:
    input:
        vcf="results/vcf/all.vcf.gz"
    output:
        vcf="results/hard_filtering/indels.vcf.gz"
    log:
        "results/logs/gatk/selectvariants/select-indels.log",
    benchmark:
        "results/benchmarks/make_indel_vcf/selectvariants-indels.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk SelectVariants -V {input.vcf}  -select-type INDEL -O {output.vcf} > {log} 2>&1 "




rule hard_filter_snps:
    input:
        vcf="results/hard_filtering/snps.vcf.gz"
    output:
        vcf="results/hard_filtering/snps-filtered.vcf.gz"
    log:
        "results/logs/gatk/variantfiltration/snps.log",
    benchmark:
        "results/benchmarks/hard_filter_snps/variantfiltration-snps.bmk"
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
        vcf="results/hard_filtering/indels.vcf.gz"
    output:
        vcf="results/hard_filtering/indels-filtered.vcf.gz"
    log:
        "results/logs/gatk/variantfiltration/indels.log",
    benchmark:
        "results/benchmarks/hard_filter_indels/variantfiltration-indels.bmk"
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
        snp="results/hard_filtering/snps-filtered.vcf.gz",
        indel="results/hard_filtering/indels-filtered.vcf.gz"
    output:
        vcf="results/vcf/all-filtered.vcf.gz",
        tbi="results/vcf/all-filtered.vcf.gz.tbi"
    log:
        "results/logs/bung_filtered_vcfs_back_together/bung.log",
    benchmark:
        "results/benchmarks/bung_filtered_vcfs_back_together/bcftools.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools concat -a {input.snp} {input.indel} | "
        " bcftools view -Oz > {output.vcf}; "
        " bcftools index -t {output.vcf}) 2> {log} "



