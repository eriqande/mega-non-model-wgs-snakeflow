

# the rules here are not ever going to be triggered
# in a regular invocation of the workflow.  This is
# for the special case of wanting to call at particular
# sites.  

# The primary use case is mapping stuff from, say, Chinook
# to Mykiss and calling it all at sites known to be SNPs
# in, say mykiss and clarkii, so that Chinook can be used
# as an outgroup for ABBA-BABA stuff.

# The basic way to do this is the include a vcf.gz file path
# as the variable `force_call_vcf` in the config YAML. A good
# place for this would be something like:
# resources/force-call/mykiss-snps.vcf.gz

# then the steps are:
# 1. Make sure that vcf.gz is indexed (.tbi)
# 2. Prepare a series of VCF files that
#   have the SNPs that we want to use for the --alleles
#   option in GATK.  This is done using the
#   scatters file.  So each of these sites.list files
#   go into force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz
#   And we index it.  
# 3. prepare a series of files that have each SNPs as a
#   region, like: chr:start-stop in a .list that GATK can
#   read, and thus limit its variant calling to those areas.
# 3. Run each of those scatters in GATK. 
# 4. bcf concat all of those.



rule index_force_call_snps:
	input:
		config["force_call_vcf"]
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/index_force_call_snps/log.stderr"
	output:
		config["force_call_vcf"] + ".tbi"
	shell:
		" bcftools index -t {input} 2> {log} "



rule prepare_alleles_option_vcfs:
	input:
		vcf=config["force_call_vcf"],
		tbi=config["force_call_vcf"] + ".tbi",
		scatters=config["scatter_intervals_file"]
	params:
		sgc="{sg_or_chrom}",
		sct="{scat}"
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/bqsr-round-{bqsr_round}/logs/prepare_alleles_option_vcfs/{sg_or_chrom}/{scat}.stderr"
	output:
		vcf="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz",
		tbi="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz.tbi"
	shell:
		" SAMP=$(bcftools query -l {input.vcf} | awk '{{print $1; exit}}' 2> {log});  "
		" REGIONS=$(awk -v id={params.sgc} -v s={params.sct} ' "
		"     $1==id && $2==s {{print $3 \":\" $4 \"-\" $5 }}  "
		"   ' {input.scatters} 2>> {log}); "
		" bcftools view -Oz -s $SAMP {input.vcf} $REGIONS > {output.vcf} 2>> {log}; "
		" bcftools index -t {output.vcf} 2>> {log}; "



rule prepare_region_files:
	input:
		vcf="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz",
		tbi="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz.tbi"
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/bqsr-round-{bqsr_round}/logs/prepare_region_files/{sg_or_chrom}/{scat}.stderr"
	output:
		reg="results/bqsr-round-{bqsr_round}/force-call/regions/{sg_or_chrom}/{scat}.list"
	shell:
		" bcftools query -f '%CHROM:%POS\n' {input.vcf} > {output.reg} 2>{log} "



rule force_call_with_gatk_scatters:
	input:
		ref="resources/genome.fasta",
		idx="resources/genome.dict",
		fai="resources/genome.fasta.fai",
		regions="results/bqsr-round-{bqsr_round}/force-call/regions/{sg_or_chrom}/{scat}.list",
		alleles="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz",
		tbi="results/bqsr-round-{bqsr_round}/force-call/alleles/{sg_or_chrom}/{scat}.vcf.gz.tbi",
		bams=expand("results/bqsr-round-{bqsr_round}/mkdup/{sample}.{ext}", bqsr_round=0, sample=sample_list, ext = ["bam"]),
		bais=expand("results/bqsr-round-{bqsr_round}/mkdup/{sample}.{ext}", bqsr_round=0, sample=sample_list, ext = ["bai"]),
	log:
		"results/bqsr-round-{bqsr_round}/logs/force_call_with_gatk_scatters/{sg_or_chrom}/{scat}.log",
	params:
		java_opts="-Xmx4g",
		hc=config["params"]["gatk"]["HaplotypeCaller"]
	conda:
		"../envs/gatk4.2.6.1.yaml"
	threads: 4
	output:
		vcf="results/bqsr-round-{bqsr_round}/force-call/vcf_sections/{sg_or_chrom}/{scat}.vcf.gz",
	resources:
		mem_mb=18800,
		time="24:00:00"
	shell:
		" if [ ! -s {input.regions} ]; then "
		" 	touch {output}; "
		" else "
		"   BAMLIST=$(echo {input.bams} | awk '{{for(i=1;i<=NF;i++) printf(\" -I %s \",$i); }}'); "
		"   gatk  --java-options \"{params.java_opts}\" HaplotypeCaller -R {input.ref} -O {output.vcf} -L {input.regions} "
		"    --alleles {input.alleles} {params.hc} --native-pair-hmm-threads {threads} "
		" $BAMLIST > {log} 2>&1; "
		" fi "


# we have to include this next rule only when the scatters file is present
if config["scatter_intervals_file"] != "":
	rule gather_gatk_force_calls:
		input:
			vcf_sections=expand("results/bqsr-round-{bqsr_round}/force-call/vcf_sections/{sg_or_chrom}/{scat}.vcf.gz",
					zip,
					bqsr_round=[0] * len(unique_scatters_table),
					sg_or_chrom=unique_scatters_table.id,
					scat=unique_scatters_table.scatter_idx,
				)
		conda:
			"../envs/bcftools.yaml"
		log:
			"results/bqsr-round-{bqsr_round}/logs/gather_gatk_force_calls/log.stderr",
		output:
			vcf="results/bqsr-round-{bqsr_round}/force-call/final.vcf.gz"
		shell:
			" (FULLS=$(for i in {input.vcf_sections}; do if [ -s $i ]; then echo $i; fi; done); "
			" bcftools concat --naive $FULLS > {output.vcf}; "
			" bcftools index -t {output.vcf} "
			" ) 2>{log} "


