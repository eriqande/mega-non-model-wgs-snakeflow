# this is the straight-up simple version that I use to just create
# a GVCF from the bam in mkdup, over all the regions.  Note that I give it
# 72 hours by default because I want it to be long enough for all possible
# bam files.
rule eca_call_variants:
    input:
        bam="results/mkdup/{sample}-{unit}.bam",
        bai="results/mkdup/{sample}-{unit}.bam.bai",
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
    output:
        gvcf=protected("results/gvcf/{sample}-{unit}.g.vcf.gz"),
        idx=protected("results/gvcf/{sample}-{unit}.g.vcf.gz.tbi"),
    conda:
        "../envs/gatk4.yaml"
    log:
        stderr="results/logs/gatk/haplotypecaller/{sample}-{unit}.stderr",
        stdout="results/logs/gatk/haplotypecaller/{sample}-{unit}.stdout",
    params:
        java_opts="-Xmx4g"
    resources:
        time="3-00:00:00",
        mem_mb = 4600,
        cpus = 1
    threads: 1
    shell:
        "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
        " -R {input.ref} "
        " -I {input.bam} "
        " -O {output.gvcf} "
        " --native-pair-hmm-threads {threads} "
        " -ERC GVCF > {log.stdout} 2> {log.stderr} "



# apparently setting Xmx4G affects the resource allocation here, but
# I want to be sure to get more RAM than that.  On slurm, it will get that
# if we request 2 CPUs on sedna.  Note that this thing runs no faster with
# four cpus and reader threads than with 2.
### NOTE: WHEN WE HAVE SAMPLE SEQUENCED IN MULTIPLE UNITS WE WILL HAVE TO COMPLETELY
### OVERHAUL THIS TO MERGE THE BAMS before Marking Duplicates, then deal just with
### the sample identifer, no longer any units, after that!!
## Note: the next two rules are unsatisfying because params.fileflags and input.gvcfs are
## defined separately, though they are effecively the same thing.
rule genomics_db_import_chromosomes:
    input:
        gvcfs=expand("results/gvcf/{u.sample}-{u.unit}.g.vcf.gz", u=units.itertuples()),
        gcvf_idxs=expand("results/gvcf/{u.sample}-{u.unit}.g.vcf.gz.tbi", u=units.itertuples()),
    output:
        db=directory("results/genomics_db/chromosomes/{chromo}"),
    log:
        "results/logs/gatk/genomicsdbimport/chromosomes/{chromo}.log"
    params:
        fileflags=expand("-V results/gvcf/{u.sample}-{u.unit}.g.vcf.gz", u=units.itertuples()),
        intervals="{chromo}",
        db_action="--genomicsdb-workspace-path", # could change to the update flag
        extra=" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --tmp-dir /scratch/eanderson/tmp ",  # optional
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenomicsDBImport {params.extra} "
        " {params.fileflags} "
        " --intervals {params.intervals} "
        " {params.db_action} {output.db} > {log} 2> {log}"



# It seems like I have to set threads to 2 in order to get 2 CPUs for this
# job.  resources: cpus = 2, does not work.  Same for above. I set reader
# threads to 2.  It doesn't seem that there is any advantage to doing more
# than that.
rule genomics_db_import_scaffold_groups:
    input:
        gvcfs=expand("results/gvcf/{u.sample}-{u.unit}.g.vcf.gz", u=units.itertuples()),
        gcvf_idxs=expand("results/gvcf/{u.sample}-{u.unit}.g.vcf.gz.tbi", u=units.itertuples()),
        scaff_groups = "scaffold_groups.tsv",
    output:
        db=directory("results/genomics_db/scaffold_groups/{scaff_group}"),
        interval_list="results/gdb_intervals/{scaff_group}.list"
    log:
        "results/logs/gatk/genomicsdbimport/scaffold_groups/{scaff_group}.log"
    params:
        fileflags=expand("-V results/gvcf/{u.sample}-{u.unit}.g.vcf.gz", u=units.itertuples()),
        db_action="--genomicsdb-workspace-path", # could change to the update flag
        extra=" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --merge-contigs-into-num-partitions 1 --tmp-dir /scratch/eanderson/tmp ",  # optional
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.yaml"
    shell:
        " export TILEDB_DISABLE_FILE_LOCKING=1; "
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output.interval_list}; "
        " gatk --java-options {params.java_opts} GenomicsDBImport {params.extra} "
        " {params.fileflags} "
        " --intervals {output.interval_list} "
        " {params.db_action} {output.db} > {log} 2> {log} "





# we can use just a single rule to create a VCF with all individuals
# for either the chromosomes or the scaffold groups. 
# {type_of_subset} will be either "chromosomes" or "scaffold_groups"
# {sg_or_chrom} will be either like "CM031199.1" (if type_of_subset is chromosomes), 
# or {scaff_group001} if type_of_subset is scaffold_groups.  
rule genomics_db2vcf:
    input:
        genome="resources/genome.fasta",
        gendb="results/genomics_db/{type_of_subset}/{sg_or_chrom}"
    output:
        vcf="results/vcf_sections/{type_of_subset}/{sg_or_chrom}.vcf.gz",
    log:
        "results/logs/gatk/genotypegvcfs/{type_of_subset}/{sg_or_chrom}.log",
    params:
        java_opts="-Xmx8g"  # I might need to consider a temp directory, too in which case, put it in the config.yaml
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "3-00:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.4.0.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        " -R {input.genome} "
        " -V gendb://{input.gendb} "
        " -O {output.vcf} > {log} 2> {log} " 


























###################################################################################
# This stuff down here is leftover from the snakemake workflow.  I had to
# do things differently, obviously.  I didn't use the wrappers in some cases
# (like GenotypeGVCFs), because that failed on the cluster in slurm mode,
# but it worked fine when I deployed it all in shell.

#    wrapper:
#        "v0.85.1/bio/gatk/genomicsdbimport"






rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        #known="resources/variation.noiupac.vcf.gz",
        #tbi="resources/variation.noiupac.vcf.gz.tbi",
        regions=(
            "called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
    wrapper:
        "0.59.0/bio/gatk/haplotypecaller"




rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand("called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index),
    output:
        gvcf="called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    wrapper:
        "0.59.2/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "0.59.2/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig=get_contigs()),
    output:
        vcf="genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    wrapper:
        "0.59.2/bio/picard/mergevcfs"
