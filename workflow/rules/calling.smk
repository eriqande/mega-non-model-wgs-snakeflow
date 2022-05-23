# this is primarily a calling thing so we will define it here:
rule make_scaff_group_interval_lists:
    input:
        scaff_groups = config["scaffold_groups"]
    output:
        "results/interval_lists/{scaff_group}.list"
    shell:
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output};"


# this is primarily a calling thing so we will define it here:
rule make_chromo_interval_lists:
    output:
        "results/interval_lists/{chromo}.list"
    shell:
        " echo {wildcards.chromo} > {output};"



# this is the straight-up simple version that I use to just create
# a GVCF from the bam in mkdup, over all the regions.  Note that I give it
# 72 hours by default because I want it to be long enough for all possible
# bam files.
# rule make_gvcfs:
#     input:
#         bam="results/mkdup/{sample}.bam",
#         bai="results/mkdup/{sample}.bai",
#         ref="resources/genome.fasta",
#         idx="resources/genome.dict",
#         fai="resources/genome.fasta.fai"
#     output:
#         gvcf=protected("results/gvcf/{sample}.g.vcf.gz"),
#         idx=protected("results/gvcf/{sample}.g.vcf.gz.tbi"),
#     conda:
#         "../envs/gatk4.2.6.1.yaml"
#     log:
#         stderr="results/logs/gatk/haplotypecaller/{sample}.stderr",
#         stdout="results/logs/gatk/haplotypecaller/{sample}.stdout",
#     benchmark:
#         "results/benchmarks/make_gvcfs/{sample}.bmk"
#     params:
#         java_opts="-Xmx4g"
#     resources:
#         time="3-00:00:00",
#         mem_mb = 4600,
#         cpus = 1
#     threads: 1
#     shell:
#         "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
#         " -R {input.ref} "
#         " -I {input.bam} "
#         " -O {output.gvcf} "
#         " --native-pair-hmm-threads {threads} "
#         " -ERC GVCF > {log.stdout} 2> {log.stderr} "



rule make_gvcf_sections:
    input:
        bam="results/mkdup/{sample}.bam",
        bai="results/mkdup/{sample}.bai",
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        fai="resources/genome.fasta.fai",
        interval_list="results/interval_lists/{sg_or_chrom}.list"
    output:
        gvcf=temp("results/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz"),
        idx=temp("results/gvcf_sections/{sample}/{sg_or_chrom}.g.vcf.gz.tbi"),
    conda:
        "../envs/gatk4.2.6.1.yaml"
    log:
        stderr="results/logs/gatk/haplotypecaller/{sample}/{sg_or_chrom}.stderr",
        stdout="results/logs/gatk/haplotypecaller/{sample}/{sg_or_chrom}.stdout",
    benchmark:
        "results/benchmarks/make_gvcfs/{sample}/{sg_or_chrom}.bmk"
    params:
        java_opts="-Xmx4g"
    resources:
        time="1-00:00:00",
        mem_mb = 4600,
        cpus = 1
    threads: 1
    shell:
        "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
        " -R {input.ref} "
        " -I {input.bam} "
        " -O {output.gvcf} "
        " -L {input.interval_list} "
        " --native-pair-hmm-threads {threads} "
        " -ERC GVCF > {log.stdout} 2> {log.stderr} "



rule concat_gvcf_sections:
    input: 
        expand("results/gvcf_sections/{{sample}}/{sgc}.g.vcf.gz", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        gvcf=protected("results/gvcf/{sample}.g.vcf.gz"),
        idx=protected("results/gvcf/{sample}.g.vcf.gz.tbi")
    log:
        "results/logs/concat_gvcf_sections/{sample}.txt"
    benchmark:
        "results/benchmarks/concat_gvcf_sections/{sample}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " bcftools concat {params.opts} -O z {input} > {output.gvcf} 2>{log}; "
        " bcftools index -t {output.gvcf} "






# apparently setting Xmx4G affects the resource allocation here, but
# I want to be sure to get more RAM than that.  On slurm, it will get that
# if we request 2 CPUs on sedna.  Note that this thing runs no faster with
# four cpus and reader threads than with 2.
## Note: the next two rules are unsatisfying because params.fileflags and input.gvcfs are
## defined separately, though they are effecively the same thing.
## This is pretty messed up.  The log file is for the last run, but I also
## tee it to one named with the import number.  
rule genomics_db_import_chromosomes:
    input:
        gvcfs=lambda wc: expand("results/gvcf/{sample}.g.vcf.gz", sample=get_samples_for_GDB_import(wc)),
        gvcf_idxs=lambda wc: expand("results/gvcf/{sample}.g.vcf.gz.tbi", sample=get_samples_for_GDB_import(wc)),
    output:
        counter="results/gdb_accounting/counters/{chromo}.txt",
        receipts=touch(expand("results/gdb_accounting/receipts/{{chromo}}/{s}", s=sample_list))
    log:
        "results/logs/gatk/genomicsdbimport/{chromo}.log"
    benchmark:
        "results/benchmarks/genomics_db_import_chromosomes/{chromo}.bmk"
    params:
        db="results/genomics_db/{chromo}",
        import_num=get_GDB_import_number,
        histories=lambda wc: expand("results/gdb_accounting/histories/{chr}/{sample}.txt", chr=wc.chromo, sample=get_samples_for_GDB_import(wc)),
        my_opts=chromo_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " mkdir -p results/genomics_db; "
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {params.db} 2>&1 | tee {log} > {log}.import_{params.import_num} && "
        " (echo $(({params.import_num} + 1)) > {output.counter}) && "
        " (for i in {params.histories}; do mkdir -p $(dirname $i); echo {params.import_num} $(date) > $i; done)"



# It seems like I have to set threads to 2 in order to get 2 CPUs for this
# job.  resources: cpus = 2, does not work.  Same for above. I set reader
# threads to 2.  It doesn't seem that there is any advantage to doing more
# than that.
rule genomics_db_import_scaffold_groups:
    input:
        gvcfs=lambda wc: expand("results/gvcf/{sample}.g.vcf.gz", sample=get_samples_for_GDB_import(wc)),
        gvcf_idxs=lambda wc: expand("results/gvcf/{sample}.g.vcf.gz.tbi", sample=get_samples_for_GDB_import(wc)),
        scaff_groups = config["scaffold_groups"],
    output:
        interval_list="results/gdb_intervals/{scaff_group}.list",
        counter="results/gdb_accounting/counters/{scaff_group}.txt",
        receipts=touch(expand("results/gdb_accounting/receipts/{{scaff_group}}/{s}", s=sample_list))
    log:
        "results/logs/gatk/genomicsdbimport/{scaff_group}.log"
    benchmark:
        "results/benchmarks/genomics_db_import_scaffold_groups/{scaff_group}.bmk"
    params:
        db="results/genomics_db/{scaff_group}",
        import_num=get_GDB_import_number,
        histories=lambda wc: expand("results/gdb_accounting/histories/{sg}/{sample}.txt", sg=wc.scaff_group, sample=get_samples_for_GDB_import(wc)),
        my_opts=scaff_group_import_gdb_opts,
        java_opts="-Xmx4g",  # optional
    resources:
        mem_mb = 9400,
        cpus = 2,
        time = "36:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " mkdir -p results/genomics_db; "
        " awk -v sg={wildcards.scaff_group} 'NR>1 && $1 == sg {{print $2}}' {input.scaff_groups} > {output.interval_list}; "
        " gatk --java-options {params.java_opts} GenomicsDBImport "
        " $(echo {input.gvcfs} | awk '{{for(i=1;i<=NF;i++) printf(\" -V %s \", $i)}}') "
        " {params.my_opts} {params.db} 2>&1 | tee {log} > {log}.import_{params.import_num} && "
        " (echo $(({params.import_num} + 1)) > {output.counter}) && "
        " (for i in {params.histories}; do mkdir -p $(dirname $i); echo {params.import_num} $(date) > $i; done) "





# {type_of_subset} will be either "chromosomes" or "scaffold_groups"
# {sg_or_chrom} will be either like "CM031199.1" (if type_of_subset is chromosomes), 
# or {scaff_group001} if type_of_subset is scaffold_groups.  
rule genomics_db2vcf:
    input:
        genome="resources/genome.fasta",
        receipts=expand("results/gdb_accounting/receipts/{{sg_or_chrom}}/{s}", s=sample_list)
    output:
        vcf="results/vcf_sections/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sections/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/gatk/genotypegvcfs/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/genomics_db2vcf/{sg_or_chrom}.bmk",
    params:
        gendb="results/genomics_db/{sg_or_chrom}",
        java_opts="-Xmx8g",  # I might need to consider a temp directory, too in which case, put it in the config.yaml
        pextra=" --genomicsdb-shared-posixfs-optimizations "
    resources:
        mem_mb = 11750,
        cpus = 2,
        time = "3-00:00:00"
    threads: 2
    conda:
        "../envs/gatk4.2.6.1.yaml"
    shell:
        " gatk --java-options {params.java_opts} GenotypeGVCFs "
        " {params.pextra} "
        " -R {input.genome} "
        " -V gendb://{params.gendb} "
        " -O {output.vcf} > {log} 2> {log} "


# this is a little rule we throw in here so that we can mark
# an individual as missing data (./. or .|.) when it has a read
# depth of 0, because GATK now marks those as 0/0,
# see https://gatk.broadinstitute.org/hc/en-us/community/posts/4476803114779-GenotypeGVCFs-Output-no-call-as-reference-genotypes?page=1#community_comment_6006727219867
# this also adds an INFO field NMISS, which gives the number of samples missing a call.
rule mark_dp0_as_missing:
    input:
        vcf="results/vcf_sections/{sg_or_chrom}.vcf.gz"
    output:
        vcf="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz",
        tbi="results/vcf_sect_miss_denoted/{sg_or_chrom}.vcf.gz.tbi"
    log:
        "results/logs/mark_dp0_as_missing/{sg_or_chrom}.log",
    benchmark:
        "results/benchmarks/mark_dp0_as_missing/{sg_or_chrom}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools +setGT {input.vcf} -- -t q -n . -i 'FMT/DP=0' | "
        " bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
        " bcftools view -Oz - > {output.vcf}; "
        " bcftools index -t {output.vcf}) 2> {log} "




rule bcf_concat:
    input:
        expand("results/hard_filtering/both-filtered-{sgc}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf=protected("results/bcf/all.bcf"),
        tbi=protected("results/bcf/all.bcf.csi")
    log:
        "results/logs/bcf_concat/bcf_concat_log.txt"
    benchmark:
        "results/benchmarks/bcf_concat/bcf_concat.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "





rule bcf_concat_mafs:
    input:
        expand("results/hard_filtering/both-filtered-{sgc}-maf-{{maf}}.bcf", sgc = unique_chromosomes + unique_scaff_groups)
    output:
        bcf=protected("results/bcf/pass-maf-{maf}.bcf"),
        tbi=protected("results/bcf/pass-maf-{maf}.bcf.csi")
    log:
        "results/logs/bcf_concat_mafs/maf-{maf}.txt"
    benchmark:
        "results/benchmarks/bcf_concat_mafs/maf-{maf}.bmk",
    params:
        opts=" --naive "
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (bcftools concat {params.opts} -Ob {input} > {output.bcf}; "
        " bcftools index {output.bcf})  2>{log}; "





