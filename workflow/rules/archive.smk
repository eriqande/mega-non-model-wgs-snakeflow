# this is a little fucntion to return the correct subd
def bam_subd(bqr):
	if bqr > 0:
		return("recal")
	else:
		return("mkdup")



# this is a simple rule that checks to make sure that the
# final requested bcfs, bams, and gvcfs have been produced.
# Then it tarballs up the qc folders and the benchmark folders
# from every bqsr-round and then it prints the rclone command that
# should be used to copy it all back to google drive.  Note that
# the BCFs from every bqsr round are copied, so that they might be
# compared later, if desired.

# this run should always be run with -np, simply to check that everything
# is done (it assumes that the qc and benchmarks are done if the other
# input is there...).  Then copy the rclone command and use it.
rule send_to_gdrive:
	input:
		expand("results/bqsr-round-{bq}/bcf/all.bcf", bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/bcf/all.bcf.csi", bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/bcf/pass-maf-{maf}.bcf", maf=mafs,  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-{fc}.txt", fc=["ALL", "PASS", "FAIL"],  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=mafs,  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/multiqc.html", bq=[str(x) for x in range(0, int(config["bqsr_rounds"])+1)]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=config["bqsr_maf"], bq=[str(x) for x in range(0, int(config["bqsr_rounds"])+1)]),
		expand("results/bqsr-round-{bq}/{subd}/{s}.{b}", bq=config["bqsr_rounds"], subd=bam_subd(config["bqsr_rounds"]), s = sample_list, b = ["bam", "bai"]),
		expand("results/bqsr-round-{bq}/gvcf/{s}.g.vcf.{ext}", bq=config["bqsr_rounds"], s=sample_list, ext=["gz", "gz.tbi"]),
		expand("results/bqsr-round-{bq}/bq_recal_tables/{s}.table", bq=[str(x) for x in range(1, int(config["bqsr_rounds"])+1)], s=sample_list),
		expand("results/bqsr-round-{bq}/bq_variants/variant{suff}", bq=[str(x) for x in range(1, int(config["bqsr_rounds"])+1)], suff=['s.vcf.gz', 's.vcf.gz.tbi', '-stats.txt']),
	params:
		BQR = config["bqsr_rounds"],
		rclone_base = config["rclone_base"],
		comma_nums = ",".join([str(x) for x in range(0,config["bqsr_rounds"]+1)]),
		bamdir="bqsr-round-{bq}/{subd}".format(bq = config["bqsr_rounds"], subd = bam_subd(config["bqsr_rounds"])),
		config_dir=os.path.dirname(config["units"]),  # assume we want to copy the whole directory that the units file is in.
		data_comm="--include='data/**'" if config["rclone_data"] else []
	shell:
		" git log | head -n 150  > {params.config_dir}/latest-git-commits.txt;  "
		" mkdir -p results/qc_summaries/bqsr-round-{{0..{params.BQR}}}; "
		" for i in {{0..{params.BQR}}}; do cp -r results/bqsr-round-$i/qc/{{multiqc.html,bcftools_stats/*.txt}} results/qc_summaries/bqsr-round-$i/; done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/qc.tar results/bqsr-round-$i/qc; gzip results/bqsr-round-$i/qc.tar;  done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/benchmarks.tar results/bqsr-round-$i/benchmarks; gzip results/bqsr-round-$i/benchmarks.tar; done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/logs.tar results/bqsr-round-$i/logs; gzip results/bqsr-round-$i/logs.tar; done; "
		"                                                                                                                  "
		" rclone copy --dry-run  --bwlimit 8650k . {params.rclone_base}  "
		" --include='{params.config_dir}/**' "
		" --include='results/qc_summaries/**' "
		" --include='results/bqsr-round-{{{params.comma_nums}}}/{{qc,benchmarks,logs}}.tar.gz' "
		" --include='results/bqsr-round-{{{params.comma_nums}}}/{{bcf,bq_recal_tables,bq_variants}}/**' "
		" --include='resources/**' "
		" {params.data_comm} "
		" --include='results/bqsr-round-{params.BQR}/gvcf/*' "
		" --include='results/{params.bamdir}/*' "
		" --include='results/bqsr-round-{params.BQR}/indel_realigned/**' "






## This version is updated to include the ANGSD-ready BAMs, etc.
## And now real BQSR (so we drop the recal stuff).
## ALSO, this won't return results from multiple rounds of BQSR
# this is a simple rule that checks to make sure that the
# final requested bcfs, bams, and gvcfs have been produced.
# Then it tarballs up the qc folders and the benchmark folders
# from every bqsr-round and then it prints the rclone command that
# should be used to copy it all back to google drive.  Note that
# the BCFs from every bqsr round are copied, so that they might be
# compared later, if desired.

# this run should always be run with -np, simply to check that everything
# is done (it assumes that the qc and benchmarks are done if the other
# input is there...).  Then copy the rclone command and use it.
rule send_to_gdrive2:
	input:
		expand("results/bqsr-round-{bq}/bcf/all.bcf", bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/bcf/all.bcf.csi", bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/bcf/pass-maf-{maf}.bcf", maf=mafs,  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-{fc}.txt", fc=["ALL", "PASS", "FAIL"],  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=mafs,  bq=config["bqsr_rounds"]),
		expand("results/bqsr-round-{bq}/qc/multiqc.html", bq=[str(x) for x in range(0, int(config["bqsr_rounds"])+1)]),
		expand("results/bqsr-round-{bq}/qc/bcftools_stats/all-pass-maf-{maf}.txt", maf=config["bqsr_maf"], bq=[str(x) for x in range(0, int(config["bqsr_rounds"])+1)]),
		expand("results/bqsr-round-{bq}/{subd}/{s}.{b}", bq=config["bqsr_rounds"], subd=bam_subd(config["bqsr_rounds"]), s = sample_list, b = ["bam", "bai"]),
		expand("results/bqsr-round-{bq}/overlap_clipped/{s}.{b}", bq=config["bqsr_rounds"], s = sample_list, b = ["bam", "bam.bai"]),
		expand("results/bqsr-round-{bq}/gvcf/{s}.g.vcf.{ext}", bq=config["bqsr_rounds"], s=sample_list, ext=["gz", "gz.tbi"]),
	params:
		BQR = config["bqsr_rounds"],
		rclone_base = config["rclone_base"],
		comma_nums = ",".join([str(x) for x in range(0,config["bqsr_rounds"]+1)]),
		bamdir="bqsr-round-{bq}/{subd}".format(bq = config["bqsr_rounds"], subd = bam_subd(config["bqsr_rounds"])),
	shell:
		" git log | head -n 150  > config/latest-git-commits.txt;  "
		" mkdir -p results/qc_summaries/bqsr-round-{{0..{params.BQR}}}; "
		" for i in {{0..{params.BQR}}}; do cp -r results/bqsr-round-$i/qc/{{multiqc.html,bcftools_stats/*.txt}} results/qc_summaries/bqsr-round-$i/; done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/qc.tar results/bqsr-round-$i/qc; gzip results/bqsr-round-$i/qc.tar;  done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/benchmarks.tar results/bqsr-round-$i/benchmarks; gzip results/bqsr-round-$i/benchmarks.tar; done; "
		" for i in {{0..{params.BQR}}}; do tar -cvf results/bqsr-round-$i/logs.tar results/bqsr-round-$i/logs; gzip results/bqsr-round-$i/logs.tar; done; "
		"                                                                                                                  "
		" rclone copy --dry-run  --drive-stop-on-upload-limit . {params.rclone_base}  "
		" --include='config/**' "
		" --include='results/qc_summaries/**' "
		" --include='results/bqsr-round-{params.BQR}/{{qc,benchmarks,logs}}.tar.gz' "
		" --include='results/bqsr-round-{params.BQR}/bcf/**' "
		" --include='resources/**' "
		" --include='data/**' "
		" --include='results/bqsr-round-{params.BQR}/gvcf/*' "
		" --include='results/{params.bamdir}/*' "
		" --include='results/bqsr-round-{params.BQR}/overlap_clipped/*' "


