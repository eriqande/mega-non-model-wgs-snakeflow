# this is the straight-up simple version that I use to just create
# a GVCF from the bam in mkdup, over all the regions.  Note that I give it
# 72 hours by default because I want it to be long enough for all possible
# bam files.
rule eca_call_variants:
    input:
        bam="results/bams_sampmerged/{sample}.bam",
        bai="results/bams_sampmerged/{sample}.bam.bai",
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        fai="resources/genome.fasta.fai"
    output:
        gvcf=protected("results/gvcf/{sample}.g.vcf.gz"),
        idx=protected("results/gvcf/{sample}.g.vcf.gz.tbi"),
    conda:
        "../envs/gatk4.2.5.0.yaml"
    log:
        stderr="results/logs/gatk/haplotypecaller/{sample}.stderr",
        stdout="results/logs/gatk/haplotypecaller/{sample}.stdout",
    benchmark:
        "results/benchmarks/haplotypecaller/{sample}.bmk"
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
        gvcfs=expand("results/gvcf/{sample}.g.vcf.gz", sample=sample_list),
        gcvf_idxs=expand("results/gvcf/{sample}.g.vcf.gz.tbi", sample=sample_list),
    output:
        db=directory("results/genomics_db/{chromo}"),
    log:
        "results/logs/gatk/genomicsdbimport/{chromo}.log"
    benchmark:
        "results/benchmarks/genomicsdbimport/{chromo}.bmk"
    params:
        fileflags=expand("-V results/gvcf/{sample}.g.vcf.gz", sample=sample_list),
        intervals="{chromo}",
        db_action="--genomicsdb-workspace-path", # could change to the update flag
        extra=" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations ",  # optional
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.5.0.yaml"
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
        gvcfs=expand("results/gvcf/{sample}.g.vcf.gz", sample=sample_list),
        gcvf_idxs=expand("results/gvcf/{sample}.g.vcf.gz.tbi", sample=sample_list),
        scaff_groups = config["scaffold_groups"],
    output:
        db=directory("results/genomics_db/{scaff_group}"),
        interval_list="results/gdb_intervals/{scaff_group}.list"
    log:
        "results/logs/gatk/genomicsdbimport/{scaff_group}.log"
    benchmark:
        "results/benchmarks/genomicsdbimport/{scaff_group}.bmk"
    params:
        fileflags=expand("-V results/gvcf/{sample}.g.vcf.gz", sample=sample_list),
        db_action="--genomicsdb-workspace-path", # could change to the update flag
        extra=" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --merge-contigs-into-num-partitions 1  ",  # optional
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.5.0.yaml"
    shell:
        " export TILEDB_DISABLE_FILE_LOCKING=1; "
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output.interval_list}; "
        " gatk --java-options {params.java_opts} GenomicsDBImport {params.extra} "
        " {params.fileflags} "
        " --intervals {output.interval_list} "
        " {params.db_action} {output.db} > {log} 2> {log} "





# {type_of_subset} will be either "chromosomes" or "scaffold_groups"
# {sg_or_chrom} will be either like "CM031199.1" (if type_of_subset is chromosomes), 
# or {scaff_group001} if type_of_subset is scaffold_groups.  
rule genomics_db2vcf:
    input:
        genome="resources/genome.fasta",
        gendb="results/genomics_db/{sg_or_chrom}"
    output:
        vcf="results/vcf_sections/{sg_or_chrom}.vcf.gz",
    log:
        "results/logs/gatk/genotypegvcfs/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/genotypegvcfs/{sg_or_chrom}.bmk",
    params:
        java_opts="-Xmx8g"  # I might need to consider a temp directory, too in which case, put it in the config.yaml
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "3-00:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.5.0.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        " -R {input.genome} "
        " -V gendb://{input.gendb} "
        " -O {output.vcf} > {log} 2> {log} " 



