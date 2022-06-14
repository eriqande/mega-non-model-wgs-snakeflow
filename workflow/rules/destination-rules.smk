# these are rules that serve the purpose of running the
# workflow to generate particular files.  They are particularly
# here so that we can easily call for the generation of files
# that need to be made so that we can check things and set parameters
# for later in the workflow.  For example, we might want to know
# the distribution of QUAL and QD values so that we can set
# a cutoff for bootstrapping the BQSR.



# I need to expand these to include all possible bqsr_rounds.  Will do that later...
rule dest_qc_0:
	input:
		expand("results/bqsr-round-0/qc/multiqc.html"),
		expand("results/bqsr-round-0/qc/bcftools_stats/all-{fc}.txt", fc=["ALL", "PASS", "FAIL"]),
		expand("results/bqsr-round-0/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=mafs),



rule dest_bqsr_histos_0:
	input:
		"results/bqsr-round-0/qc/bqsr_relevant_histograms/qd.tsv",
		"results/bqsr-round-0/qc/bqsr_relevant_histograms/qual.tsv"