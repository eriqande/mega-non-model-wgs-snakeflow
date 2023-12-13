import os
import warnings
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

# The two lines below can be used to get the config information
# for testing.  Change config file path as appropriate
#  from snakemake.io import _load_configfile
#  config = _load_configfile("example-configs/nz-arg-wr-chinook/config.yaml")

min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
#configfile: config["config"]



validate(config, schema="../schemas/config.schema.yaml")


# these are our MAF cutoffs to prepare.  We take the unique values
# of the union of the maf_cutoffs from the config and the bqsr_maf
# from the config.
mafs = list(
        dict.fromkeys(
            [str(x) for x in config["maf_cutoffs"]] + 
            [str(config["bqsr_maf"])]
            )
        )



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

# this is handy for getting sample info from the bcftools summaries
unique_sample_ids = list(units["sample_id"].unique())

### Eric's addition to genotype over the chromosomes and scaffold_groups
chromosomes = pd.read_table(config["chromosomes"]).set_index("chrom", drop=False)
validate(chromosomes, schema="../schemas/chromosomes.schema.yaml")

scaffold_groups = pd.read_table(config["scaffold_groups"]).set_index("id", drop=False)
validate(scaffold_groups, schema="../schemas/scaffold_groups.schema.yaml")

# ensure that column 1 of the scaffold group file is "id" and
# column 2 is "chrom".  This is essential because I used those
# column positions in some awk to pull things out.
scaff_cols = list(scaffold_groups.columns)
if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom': 
    raise Exception("Column order is important in the scaffold_groups file.  The first column must be 'id' and the second column must be 'chrom'.")

# get a list of just the unique values of the scaffold_group and of the chromosomes
unique_scaff_groups = list(scaffold_groups.id.unique())
unique_chromosomes = list(chromosomes.chrom.unique())  # don't need to unique it, but I do anyway


# finally, get all the scatter groups, indexed two different ways.
scatter_wc_constraint="scat_0[0-9]*"  # this is just here for the case where `scatter_intervals_file: ""`
if config["scatter_intervals_file"] != "":
    scatter_groups = pd.read_table(config["scatter_intervals_file"]).set_index("id", drop=False)
    validate(scatter_groups, schema="../schemas/scatter_intervals.schema.yaml")
    scatter_cols = list(scatter_groups.columns)
    if scatter_cols[0] != 'id' or scatter_cols[1] != 'scatter_idx' or scatter_cols[2] != 'chrom' or scatter_cols[3] != 'start' or scatter_cols[4] != 'end' or scatter_cols[5] != 'scatter_length':
        raise Exception("Column order is important in the scaffold_groups file.  The columns must be in order: id, scatter_idx, chrom, start, end, scatter_length.")
    unique_scats = list(scatter_groups.scatter_idx.unique())
    scatter_wc_constraint="|".join(unique_scats)
    # get a pandas frame of unique values of scaff_group and scat. This is for force-calling VCFs
    unique_scatters_table=scatter_groups[['id', 'scatter_idx']].drop_duplicates()





### Deal with the indel_grps if present (i.e. groupings of the samples
### into different species so that indel realignment is done species by species).
### At the same time we do this, we are also going to define the lists of output bams

# this is the default, but we update it if the indel_grps is defined in the config
realigned_bams_output_list=expand("results/bqsr-round-{bq}/indel_realigned/__ALL/{samp}.bam", bq=config["bqsr_rounds"], samp=sample_list)
indel_grps_list=["weirdwhacko", "__ALL"]  # this is for the wildcard constraints. For some reason it fails if it is just "__ALL"
if "indel_grps" in config and config["indel_grps"] != "":
    indel_grps = pd.read_table(config["indel_grps"], dtype=str).set_index("sample", drop=False)
    validate(indel_grps, schema="../schemas/indel_grps.schema.yaml")

    # check to make sure that every sample is accounted for
    if(set(sample_list) != set(indel_grps["sample"].unique())):
      raise Exception("Every sample in units must be represented in the indel_grps file")

    # update the realigned_bams_output_list
    sm=indel_grps['sample'].tolist()
    ig=indel_grps['indel_grp'].tolist()
    realigned_bams_output_list=["results/bqsr-round-{bq}/indel_realigned/{ig}/{samp}.bam".format(bq=config["bqsr_rounds"], ig=ig[i], samp=sm[i]) for i in range(len(ig))]
    indel_grps_list=ig


# down here, if we choose not to do indel realignment in the config
# then we just set realigned_bams_output_list to an empty list
if "do_indel_realignment" in config and config["do_indel_realignment"] == False:
    realigned_bams_output_list = []


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(sample_list),
    unit="|".join(units["unit"]),
    chromo="|".join(unique_chromosomes),
    scaff_group="|".join(unique_scaff_groups),
    sg_or_chrom="|".join(unique_scaff_groups + unique_chromosomes),
    filter_condition="ALL|PASS|FAIL",
    maf="|".join(mafs),
    scatter=scatter_wc_constraint,
    igrp="|".join(indel_grps_list),
    bqsr_round="|".join(["0","1","2","3","4"])



#### Pick out all the units that are of the same sample in the same library
# def get_units_of_common_sample_and_lib(wildcards):
#     su = units.loc[(units["sample"] == wildcards.sample) & (units["library"] == wildcards.library)]
#     return(expand("results/bqsr-round-{bqsr_round}/mapped/{sample}---{unit}.sorted.bam", zip,
#         sample = su["sample"].tolist(),
#         unit = su["unit"].tolist(),
#     ))

# get all the units of a particular sample
def get_all_bams_of_common_sample(wildcards):
    s=units.loc[(units["sample"] == wildcards.sample)]
    return(expand("results/bqsr-round-{{bqsr_round}}/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist(),
    ))


# #### Pick out all the libmerged files of a sample
# def get_libmerged_bams_of_common_sample(wildcards):
#     su = units.loc[(units["sample"] == wildcards.sample)]
#     # make a list of all libmerged bams
#     dupie_list = expand(
#         "results/bqsr-round-{bqsr_round}/mkdup/{sample}---{library}.bam",
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
# results/bqsr-round-{bqsr_round}/gdb_accounting/receipts/{chrom_or_scaff_group}/sample_name.
# The contents of the file are a number
# which is the number of the update.  0 = initial import, 1 = first
# update, etc.  and on the same line we will put the date that happened.
#
# Additionally, each time the data base is imported or updated, we will
# write a line with the total number of imports/updates made to
# results/bqsr-round-{bqsr_round}/gdb_accounting/counters/{chrom_or_scaff_group.txt}.
# We add those up to get what the number should be.

# here chrom is either a chromosome or a scaffold group, because we have to get
# either of those.
def previously_imported_samples(chrom, bq):
    dir="results/bqsr-round-{bq}/gdb_accounting/receipts/".format(bq = bq) + chrom
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
    return(list(set(sample_list).difference(set(previously_imported_samples(chrom, wildcards.bqsr_round)))))


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
    file="results/bqsr-round-{{bqsr_round}}/gdb_accounting/counters/" + chrom + ".txt"
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
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --intervals results/bqsr-round-{bq}/gdb_intervals/{sg}.list --merge-contigs-into-num-partitions 1  --genomicsdb-workspace-path ".format(bq = wildcards.bqsr_round, sg = wildcards.scaff_group))
    else:
        return(" --batch-size 50 --reader-threads 2 --genomicsdb-shared-posixfs-optimizations --merge-contigs-into-num-partitions 1  --genomicsdb-update-workspace-path ")



###################################################################################################

## Here we get the -L option(s) for a chromosome or a scaff_group

## Here is one for figuring out how to filter for bcftools stats
def get_bcftools_stats_filter_option(wildcards):
    if(wildcards.filter_condition == "ALL"):
        return(" ")
    elif(wildcards.filter_condition == "PASS"):
        return( " -i 'FILTER=\"PASS\"' ")
    elif(wildcards.filter_condition == "FAIL"):
        return( " -i 'FILTER!=\"PASS\"' ")
    else:
        raise Exception("Wildcard filter_condition must be ALL, PASS, or FAIL.")




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


# here is the function that picks out the appropriate bam and bai file
# to use for calling gvcfs.  If the {bqsr_round} is = 0 it pulls them
# from bqsr_round-{bqsr-round}/mkdup. 
# If the {bqsr_round} is > 0 it pulls them
# from bqsr_round-{bqsr-round}/recal.
def get_bams_for_calling(wildcards):
    if wildcards.bqsr_round == "0":
        subd = "mkdup"
    else:
        subd = "recal"
    return { 
        "bam": "results/bqsr-round-{bqsr_round}/{subd}/{sample}.bam".format(
            bqsr_round = wildcards.bqsr_round,
            subd = subd,
            sample = wildcards.sample),
        "bai": "results/bqsr-round-{bqsr_round}/{subd}/{sample}.bai".format(
            bqsr_round = wildcards.bqsr_round,
            subd = subd,
            sample = wildcards.sample)}

def get_bams_for_samtools_stats(wildcards):
    if wildcards.bqsr_round == "0":
        subd = "mkdup"
    else:
        subd = "recal"
    return "results/bqsr-round-{bqsr_round}/{subd}/{sample}.bam".format(
            bqsr_round = wildcards.bqsr_round,
            subd = subd,
            sample = wildcards.sample)

# this function looks back to the previous round of BQSR and
# uses the resulting MAF file from there as input
def get_variants_to_condense_for_bqsr(wildcards):
    bqr = int(wildcards.bqsr_round) - 1
    return("results/bqsr-round-{bq}/bcf/pass-maf-{maf}.bcf".format(bq=bqr, maf = config["bqsr_maf"]))


def get_bams_for_bqsr(wildcards):
    if wildcards.bqsr_round == "0":
        return("");
    elif wildcards.bqsr_round == "1":
        subd="mkdup"
    else: 
        subd="recal"
    bqr = int(wildcards.bqsr_round) - 1
    return { 
        "bam": "results/bqsr-round-{bqsr_round}/{subd}/{sample}.bam".format(
            bqsr_round = bqr,
            subd = subd,
            sample = wildcards.sample),
        "bai": "results/bqsr-round-{bqsr_round}/{subd}/{sample}.bai".format(
            bqsr_round = bqr,
            subd = subd,
            sample = wildcards.sample)}


# given a chromosome or a scaff_group, get all the scattered vcf.gz of vcf.gz.tbi files
def get_scattered_vcfs(wildcards, ext):
    scat_ids=scatter_groups.loc[(scatter_groups["id"] == wildcards.sg_or_chrom), "scatter_idx"].unique()
    return expand("results/bqsr-round-{{bqsr_round}}/vcf_sections/{{sg_or_chrom}}/{scat}.vcf.gz{e}", scat=scat_ids, e=ext)



# we have this here becuase we only want to do fastqc, mkdup and trimmomatic
# qc for bqsr_round=0.  The others just do the samtools stats.
def get_multiqc_inputs(wildcards):
    if wildcards.bqsr_round == "0":
        return(
            expand("results/bqsr-round-{bq}/qc/samtools_stats/{sample}.txt", bq=wildcards.bqsr_round, sample=sample_list) +
            expand("results/bqsr-round-{bq}/qc/fastqc/{u.sample}---{u.unit}_R{r}_fastqc.zip", bq=wildcards.bqsr_round, u=units.itertuples(), r = [1, 2]) +
            list(dict.fromkeys(expand("results/bqsr-round-{bq}/qc/mkdup/{u.sample}.metrics.txt", bq=wildcards.bqsr_round, u=units.itertuples()))) +
            expand("results/bqsr-round-{bq}/logs/trim_reads_pe/{u.sample}---{u.unit}.log", bq=wildcards.bqsr_round, u=units.itertuples())
        )
    else:
        return expand("results/bqsr-round-{bq}/qc/samtools_stats/{sample}.txt", bq=wildcards.bqsr_round, sample=sample_list)
    

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/bqsr-round-{{bqsr_round}}/trimmed/{sample}---{unit}.{group}.fastq.gz", group=[1, 2], **wildcards
        )
    # single end sample
    return "results/bqsr-round-{{bqsr_round}}/trimmed/{sample}---{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/bqsr-round-{{bqsr_round}}/mkdup/{sample}---{unit}.bam",
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
    f = "results/bqsr-round-{bqsr_round}/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/bqsr-round-{bqsr_round}/mkdup/{sample}-{unit}.bam"
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


# here, we deal with indel realignment by species (or "igrp")
def get_igrp_sample_names(wildcards):
    if "indel_grps" not in config or config["indel_grps"] == "":
        if wildcards.igrp != "__ALL":
            raise Exception("Requesting igrp that is not \"__ALL\", but indel_grps not defined in config.")
        else:
            ret="JustStuffNotNeeded: IDs gotten through Unix in the rule if __ALL"
    else:
        if wildcards.igrp == "__ALL":
          raise Exception("You can't request igrp \"__ALL\" when indel_grps is defined in the config. ")
        else:
          ret=indel_grps.loc[indel_grps['indel_grp'] == wildcards.igrp]["sample_id"].unique().tolist()
    return ret


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    )


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}
