


# make a VCF file that only has one of the samples in it
# for GATK to use
rule condense_variants_for_bqsr:
	input: 
		bcf = get_variants_to_condense_for_bqsr,
		units = config["units"]
	output: 
		vcf = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz",
		tbi = "results/bqsr-round-{bqsr_round}/bq_variants/variants.vcf.gz.tbi"
	log:
		"results/bqsr-round-{bqsr_round}/logs/condense_variants_for_bqsr/bcftools.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"SAMP=$(bcftools query -l all.bcf | head -n 1); "
		" bcftools view -s $SAMP {input.bcf} > {output.vcf} 2> {log}; "
		" bcftools index -t {output} 2>> {log}; "
