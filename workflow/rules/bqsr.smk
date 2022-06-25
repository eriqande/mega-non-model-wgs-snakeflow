


# This is a rule that will create histogram summaries of the
# distribution of QUALS and CDs.  These can be used to
# set the bqsr values.
rule bqsr_relevant_histograms:
	input:
		bcf = "results/bqsr-round-{{bqsr_round}}/bcf/pass-maf-{maf}.bcf".format(maf = config["bqsr_maf"])
	output:
		qd = "results/bqsr-round-{bqsr_round}/qc/bqsr_relevant_histograms/qd.tsv",
		qual = "results/bqsr-round-{bqsr_round}/qc/bqsr_relevant_histograms/qual.tsv"
	log:
		"results/bqsr-round-{bqsr_round}/logs/bqsr_relevant_histograms/log.txt"
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/bqsr_relevant_histograms/benchmark.bmk"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"workflow/scripts/qd_and_qual.sh {input.bcf} {output.qd} {output.qual} {log} "


# make a VCF file that only has one of the samples in it
# for GATK to use
rule condense_variants_for_bqsr:
	input: 
		bcf = get_variants_to_condense_for_bqsr,
		units = config["units"]
	output: 
		vcf = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz",
		tbi = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz.tbi",
		stats = "results/bqsr-round-{bqsr_round}/bq_variants/variant-stats.txt",
	params:
		qual = config["bqsr_qual"],
		qd = config["bqsr_qd"]
	log:
		"results/bqsr-round-{bqsr_round}/logs/condense_variants_for_bqsr/bcftools.log"
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/condense_variants_for_bqsr/benchmark.bmk"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"SAMP=$(bcftools query -l {input.bcf} | head -n 1); "
		" bcftools view -i 'QUAL >= {params.qual} && QD >= {params.qd}' "
		"     -Oz -s $SAMP {input.bcf} > {output.vcf} 2> {log}; "
		" bcftools index -t {output.vcf} 2>> {log}; "
		" bcftools stats {output.vcf} > {output.stats} 2>> {log}; "



rule recalibrate_bases:
	input: 
		unpack(get_bams_for_bqsr),  # this gets us the bam and bai for the sample
		vcf = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz",
		tbi = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz.tbi",
		stats = "results/bqsr-round-{bqsr_round}/bq_variants/variant-stats.txt", # it doesn't need the stats, but this will at least trigger their creation
		ref="resources/genome.fasta",
	output:
		"results/bqsr-round-{bqsr_round}/bq_recal_tables/{sample}.table"
	log:
		"results/bqsr-round-{bqsr_round}/logs/gatk/baserecalibrator/{sample}.log"
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/recalibrate_bases/{sample}.bmk"
	conda:
		"../envs/gatk4.2.6.1.yaml"
	shell:
		"gatk BaseRecalibrator "
		" -I {input.bam} "
		" -R {input.ref} "
		" --known-sites {input.vcf} "
		" --bqsr-baq-gap-open-penalty 30 "
		" -O {output}  2> {log} "




rule apply_bqsr:
	input: 
		unpack(get_bams_for_bqsr),  # this gets us the bam and bai for the sample
		ref="resources/genome.fasta",
		rtable="results/bqsr-round-{bqsr_round}/bq_recal_tables/{sample}.table"
	output:
		bam="results/bqsr-round-{bqsr_round}/recal/{sample}.bam",
		bai="results/bqsr-round-{bqsr_round}/recal/{sample}.bai"
	log:
		"results/bqsr-round-{bqsr_round}/logs/gatk/applybqsr/{sample}.log"
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/apply_bqsr/{sample}.bmk"
	conda:
		"../envs/gatk4.2.6.1.yaml"
	shell:
		"gatk --java-options \"-Dsamjdk.compression_level=9\" ApplyBQSR "
		" -R {input.ref} "
		" -I {input.bam} "
		" --bqsr-recal-file {input.rtable} "
		" -O {output.bam}  2> {log} "
