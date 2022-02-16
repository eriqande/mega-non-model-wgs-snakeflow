import glob
import re
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version



min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


### Eric's addition to genotype over the chromosomes and scaffold_groups
chromosomes = pd.read_table("chromosomes.tsv").set_index("id", drop=False)
validate(chromosomes, schema="../schemas/chromosomes.schema.yaml")

scaffold_groups = pd.read_table("scaffold_groups.tsv").set_index("id", drop=False)
validate(scaffold_groups, schema="../schemas/scaffold_groups.schema.yaml")

# get a list of just the unique values of the scaffold_group
unique_scaff_groups = list(scaffold_groups.id.unique())


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Eric's Helper Functions to Handle GenomicDB updating

# here chrom is either a chromosome or a scaffold group, because we have to get
# either of those.
def check_dbi_receipts(chrom):
    # get the different files in two lists
    adds_receipts = glob.glob("results/dbi_run_receipts/" + chrom + "/*/added_samples.txt", recursive=True)
    cums_receipts = glob.glob("results/dbi_run_receipts/" + chrom + "/*/cumulative_samples.txt", recursive=True)
    # extract the integer part
    adds_nums = [int(re.findall("[0-9]+", i )[0])  for i in adds_receipts]
    cums_nums = [int(re.findall("[0-9]+", i )[0])  for i in cums_receipts]
    adds_nums.sort()
    cums_nums.sort()
    # make sure that the two files are present in all numbered subdirs
    if(adds_nums != cums_nums):
        raise Exception("added_samples.txt and cumulative_samples.txt occurrences not congruent in results/dbi_run_receipts/[0-9]+/")
    # if there are no such files, we return a -1
    if not bool(adds_receipts):
        return(-1)
    # get the maximum value in adds_nums are return it
    return(max(adds_nums))


def get_samples_for_GDB_import(wildcards):
    # check for receipts and corrctness thereof
    lastDBImport = check_dbi_receipts()
    # get the written number of the run from the config
    gdi_run = config["genomics_db_import_num"]
    if(gdi_run > lastDBImport + 1):
        raise Exception("Requested value of config variable `genomics_db_import_num` must be, at most, one larger than the highest dbi_run_receipt.")


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
    return r"-R '@RG\tID:{sample}\tSM:{sample_id}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}'".format(
        sample=wildcards.sample,
        sample_id=units.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=units.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=units.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=units.loc[(wildcards.sample, wildcards.unit), "lane"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/mkdup/{sample}-{unit}.bam",
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
