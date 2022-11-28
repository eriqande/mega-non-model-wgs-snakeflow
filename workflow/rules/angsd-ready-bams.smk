
# These are rules that take the mapped BAMs, and the results
# of running GATK to get the indels, and do:
#    1. Clipping over overlaps
#    2.  Indel realignment
# 
# This is set up to allow these things to happen for any
# bqsr_round, but the truth is, we will probably simply use
# it as bqsr_round 0, unless we start to get crazy-complete
# variant data for our non-model organisms...highly unlikely




# Simple clipping of overlaps in the mkduped BAMs
rule clip_overlaps:
	input:
		"results/bqsr-round-{bqsr_round}/mkdup/{sample}.bam"
	output:
		bam="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam",
		bai="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam.bai"
	log:
		clip="results/bqsr-round-{bqsr_round}/logs/clip_overlaps/clip_overlap-{sample}.log",
		index="results/bqsr-round-{bqsr_round}/logs/clip_overlaps/index-{sample}.log"
	conda:
		"../envs/bamutil_samtools.yaml"
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/clip_overlaps/{sample}.bmk"
	shell:
		" bam clipOverlap --in {input} --out {output} --stats 2> {log.clip} && "
		" samtools index {output} 2> {log.index}"


rule species_sample_lists:
    params:
        sams=get_igrp_sample_names
    output:
        sample_list="results/igrp_lists/{igrp}.args",
    log:
        "results/logs/species_sample_lists/{igrp}.log"
    shell:
        " for i in {params.sams}; do echo $i; done > {output.sample_list} 2>{log} "




# When we have multiple species that we are dealing with,
# then we have a step of making species specific indel VCFs
rule make_species_specific_indel_vcfs:
    input:
        indel="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz",
        indel_idx="results/bqsr-round-{bqsr_round}/hard_filtering/indels-filtered-{sg_or_chrom}.vcf.gz.tbi",
        sample_list="results/igrp_lists/{igrp}.args",
        ref="resources/genome.fasta"
    output:
         vcf=temp("results/bqsr-round-{bqsr_round}/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz"),
         idx=temp("results/bqsr-round-{bqsr_round}/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz.tbi"),
    log:
        "results/bqsr-round-{bqsr_round}/logs/make_species_specific_indel_vcfs/{igrp}/{sg_or_chrom}.log",
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/make_species_specific_indel_vcfs/{igrp}/selectvariants-{sg_or_chrom}.bmk"
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " IGRP={wildcards.igrp};                           "  
        " if [ $IGRP = \"__ALL\" ]; then                   "   # just hard-link the files in this case
        "    (ln  {input.indel} {output.vcf};              "
        "    ln  {input.indel_idx} {output.idx}) 2>{log};  "
        " else                                        "
        "   gatk SelectVariants -R {input.ref}        "
        "     -V {input.indel}                        "
        "     -O {output.vcf}                         "
        "     --sample-name {input.sample_list}       "
        "     --exclude-non-variants 2>{log};         "
        " fi                                          "

# We need to get a VCF file holding all the known indels.
# We will use the ones that we hardfiltered.  We only keep
# the first sample in each to make it a lot smaller.  (We
# could probably make a VCF with no samples in it, but it
# is easy enough to just keep one sample, and we know that
# will work).  We do this on a chromosome or scaff group basis
rule get_known_indels:
	input:
		indel="results/bqsr-round-{bqsr_round}/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz",
		indel_idx="results/bqsr-round-{bqsr_round}/species-spec-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz.tbi"
	output:
		vcf=temp("results/bqsr-round-{bqsr_round}/known-indel-sections/{igrp}/{sg_or_chrom}.vcf.gz"),
	log:
		"results/bqsr-round-{bqsr_round}/logs/get_known_indels/{igrp}/{sg_or_chrom}.log",
	benchmark:
		"results/bqsr-round-{bqsr_round}/benchmarks/get_known_indels/bcftools-{igrp}-{sg_or_chrom}.bmk"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"SAMP=$(bcftools query -l {input.indel} | head -n 1); "
		" bcftools view -Oz -s $SAMP {input.indel} > {output.vcf} 2> {log}; "
		" bcftools index -t {output.vcf} 2>> {log}; "


# now we concat the above outputs all together into a singl vcf
rule bcfconcat_known_indels:
    input:
    	expand("results/bqsr-round-{{bqsr_round}}/known-indel-sections/{{igrp}}/{sgc}.vcf.gz", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        vcf=protected("results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz"),
        idx=protected("results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz.tbi")
    log:
        "results/bqsr-round-{bqsr_round}/logs/bcfconcat_known_indels/{igrp}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/bcfconcat_known_indels/bcf_concat-{igrp}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Oz {input} > {output.vcf} 2> {log}; "
        " bcftools index -t {output.vcf})  2>> {log}; "


# Download the gatk3.8 jar and then run gatk3-register. This is necessary
# because the licensing for gatk3 doesn't allow the jar to be distributed
# directly off of conda.  Apparently this only has to be done for
# Mac OS X---gatk3-register is not even available on the Linux
# conda package.  What-evs.
rule gatk3_register:
    output:
    	dir=directory("results/gatk3.8_jar"),
    	jar="results/gatk3.8_jar/GenomeAnalysisTK.jar"
    log:
        "results/logs/gatk3_register/log.txt"
    benchmark:
        "results/benchmarks/gatk3_register/gatk3_register.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
    	jopts="-Xmx4g"
    shell:
        "( wget -P {output.dir} https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 && "
        "  bzip2 -d  {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2  && "
        "  tar --directory {output.dir} -xvf {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836.tar && "
        "  mv {output.dir}/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar {output.dir} && "
        "  if [ $(uname -s) = \"Darwin\" ]; then gatk3-register {output.jar}; fi"
        "  ) > {log} 2>&1 "





# Run RealignerTargetCreator from GATK 3.8.0.
rule realigner_target_creator:
    input:
    	vcf="results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz",
    	idx="results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz.tbi",
    	jar="results/gatk3.8_jar/GenomeAnalysisTK.jar"
    output:
    	"results/bqsr-round-{bqsr_round}/realigner-targets/{igrp}/realigner-targets.intervals"   
    log:
        "results/bqsr-round-{bqsr_round}/logs/realigner_target_creator/{igrp}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/realigner_target_creator/{igrp}.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
    	jopts="-Xmx4g"
    threads: 4
    shell:
        "gatk3 {params.jopts} -T RealignerTargetCreator -nt {threads} "
        " -R resources/genome.fasta "
        " -o {output} "
        " -known {input.vcf} 2> {log} "



# realign indels on each bam using the known indel variation
# stored in the realigner_target_creator step
rule indel_realigner:
    input:
    	vcf="results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz",
    	idx="results/bqsr-round-{bqsr_round}/known-indels/{igrp}/all-indels.vcf.gz.tbi",
    	jar="results/gatk3.8_jar/GenomeAnalysisTK.jar",
    	intervals="results/bqsr-round-{bqsr_round}/realigner-targets/{igrp}/realigner-targets.intervals",
    	fasta="resources/genome.fasta",
    	bam="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam",
    	bai="results/bqsr-round-{bqsr_round}/overlap_clipped/{sample}.bam.bai"
    output:
    	bam="results/bqsr-round-{bqsr_round}/indel_realigned/{sample}-{igrp}.bam"
    log:
        "results/bqsr-round-{bqsr_round}/logs/indel_realigner/{sample}-{igrp}.txt"
    benchmark:
        "results/bqsr-round-{bqsr_round}/benchmarks/indel_realigner/{sample}-{igrp}.bmk",
    conda:
        "../envs/gatk3.8.yaml"
    params:
    	jopts="-Xmx4g"
    shell:
        "gatk3 {params.jopts}  \"-Dsamjdk.compression_level=9\" "
        " -T IndelRealigner  "
        " -R {input.fasta} "
        " -I {input.bam} "
        " -targetIntervals {input.intervals} "
        " -known {input.vcf} "
        " --consensusDeterminationModel KNOWNS_ONLY "
        " -LOD 0.4 " # this is recommended at https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Local_Realignment_around_Indels.md
        " -o {output.bam} "
        " > {log} 2>&1 "





