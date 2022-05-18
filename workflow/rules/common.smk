import os
import warnings
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version



min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


# rather than have a separate samples.tsv, we can just get a list of
# the samples from units
sample_list = list(units["sample"].unique())

### Eric's addition to genotype over the chromosomes and scaffold_groups
chromosomes = pd.read_table(config["chromosomes"]).set_index("chrom", drop=False)
validate(chromosomes, schema="../schemas/chromosomes.schema.yaml")

scaffold_groups = pd.read_table(config["scaffold_groups"]).set_index("id", drop=False)
validate(scaffold_groups, schema="../schemas/scaffold_groups.schema.yaml")

# get a list of just the unique values of the scaffold_group and of the chromosomes
unique_scaff_groups = list(scaffold_groups.id.unique())
unique_chromosomes = list(chromosomes.chrom.unique())  # don't need to unique it, but I do anyway


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(sample_list),
    unit="|".join(units["unit"]),
    chromo="|".join(chromosomes["chrom"]),
    scaff_group="|".join(unique_scaff_groups),
    sg_or_chrom="|".join(unique_scaff_groups + list(chromosomes["chrom"]))


#### Pick out all the units that are of the same sample in the same library
# def get_units_of_common_sample_and_lib(wildcards):
#     su = units.loc[(units["sample"] == wildcards.sample) & (units["library"] == wildcards.library)]
#     return(expand("results/mapped/{sample}---{unit}.sorted.bam", zip,
#         sample = su["sample"].tolist(),
#         unit = su["unit"].tolist(),
#     ))

# get all the units of a particular sample
def get_all_bams_of_common_sample(wildcards):
    s=units.loc[(units["sample"] == wildcards.sample)]
    return(expand("results/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist(),
    ))


# #### Pick out all the libmerged files of a sample
# def get_libmerged_bams_of_common_sample(wildcards):
#     su = units.loc[(units["sample"] == wildcards.sample)]
#     # make a list of all libmerged bams
#     dupie_list = expand(
#         "results/mkdup/{sample}---{library}.bam",
#         zip,
#         sample = su["sample"].tolist(),
#         library = su["library"].tolist(),
#     )
#     # then return just the unique elements from that
#     return(list(dict.fromkeys(dupie_list)))



##################################################################
##### Eric's Helper Functions to Handle GenomicDB updating  #####

# here is how we are doing it.  Every time a sample gets put
# in the data base, a file with the sample name is put into
# results/gdb_accounting/receipts/{chrom_or_scaff_group}/sample_name.
# The contents of the file are a number
# which is the number of the update.  0 = initial import, 1 = first
# update, etc.  and on the same line we will put the date that happened.
#
# Additionally, each time the data base is imported or updated, we will
# write a line with the total number of imports/updates made to
# results/gdb_accounting/counters/{chrom_or_scaff_group.txt}.
# We add those up to get what the number should be.

# here chrom is either a chromosome or a scaffold group, because we have to get
# either of those.
def previously_imported_samples(chrom):
    dir="results/gdb_accounting/receipts/" + chrom
    if not  os.path.exists(dir):
        return([])  
    return(os.listdir(dir))
    


# samples to put in the GDB are any in sample list
# that have not already been imported.
def get_samples_for_GDB_import(wildcards):
    # wildcards might contain either scaff_group or chromo.
    # whichever it is, we get it and save it as the variable chrom
    if hasattr(wildcards, 'chromo'):
        chrom=wildcards.chromo
    elif hasattr(wildcards, 'scaff_group'):
            chrom=wildcards.scaff_group
    else:
        raise Exception("Wildcards has neither chromo or scaff_group in get_samples_for_GDB_import")
    return(list(set(sample_list).difference(set(previously_imported_samples(chrom)))))


# get the number of data base import or update
# based on the dbi_counters
def get_GDB_import_number(wildcards):
    # wildcards might contain either scaff_group or chromo.
    # whichever it is, we get it and save it as the variable chrom
    if hasattr(wildcards, 'chromo'):
        chrom=wildcards.chromo
    elif hasattr(wildcards, 'scaff_group'):
            chrom=wildcards.scaff_group
    else:
        raise Exception("Wildcards has neither chromo or scaff_group in get_samples_for_GDB_import")
    file="results/gdb_accounting/counters/" + chrom + ".txt"
    if not  os.path.exists(file):
        return(0)  
    return int(open(file).readlines()[0])




def chromo_import_gdb_opts(wildcards):
    inum=get_GDB_import_number(wildcards)
    if inum == 0:
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --intervals {chr} --genomicsdb-workspace-path ".format(chr = wildcards.chromo))
    else:
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --genomicsdb-update-workspace-path ")


def scaff_group_import_gdb_opts(wildcards):
    inum=get_GDB_import_number(wildcards)
    if inum == 0:
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --intervals results/gdb_intervals/{sg}.list --merge-contigs-into-num-partitions 1  --genomicsdb-workspace-path ".format(sg = wildcards.scaff_group))
    else:
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --merge-contigs-into-num-partitions 1  --genomicsdb-update-workspace-path ")


###################################################################################################



##### Helper functions #####
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{sample_id}_{library}_{flowcell}_{lane}_{barcode}\tSM:{sample_id}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}.{barcode}'".format(
        sample=wildcards.sample,
        sample_id=units.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=units.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=units.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=units.loc[(wildcards.sample, wildcards.unit), "lane"],
        barcode=units.loc[(wildcards.sample, wildcards.unit), "barcode"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}---{unit}.{group}.fastq.gz", group=[1, 2], **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}---{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/mkdup/{sample}---{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/mkdup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    )


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}
