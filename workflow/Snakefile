include: "rules/common.smk"


##### Target rules #####


rule all:
    input:
        #expand("results/genomics_db/chromosomes/{c}", c=chromosomes.chrom),
        #expand("results/genomics_db/scaffold_groups/{s}", s=unique_scaff_groups),
        expand("results/vcf_sections/chromosomes/{c}.vcf.gz", c=chromosomes.chrom),
        expand("results/vcf_sections/scaffold_groups/{s}.vcf.gz", s=unique_scaff_groups),
        "results/qc/multiqc.html",

rule old_all:
    input:
        "annotated/all.vcf.gz",
        "qc/multiqc.html",
        "plots/depths.svg",
        "plots/allele-freqs.svg",


##### Modules #####


include: "rules/ref.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
