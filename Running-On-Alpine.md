Running this Workflow on Alpine
================
Eric C. Anderson

- [Fresh install of mambaforge](#fresh-install-of-mambaforge)
- [Get Snakemake](#get-snakemake)
- [Get mega-non-model-wgs snakeflow and test it in an interactive
  session](#get-mega-non-model-wgs-snakeflow-and-test-it-in-an-interactive-session)
- [Start copying data over](#start-copying-data-over)
- [copy the old resources directory that has the indexed
  genome](#copy-the-old-resources-directory-that-has-the-indexed-genome)
- [Make an alpine slurm profile](#make-an-alpine-slurm-profile)
- [Running](#running)
- [Restarting](#restarting)
- [Reboot! The sample IDs were not the same between
  runs](#reboot-the-sample-ids-were-not-the-same-between-runs)

This documents Eric modifying the workflow to run on alpine, and also
getting things setup to run all the SWTH data.

# Fresh install of mambaforge

To run this workflow on Alpine, I wanted to start from a completely
clean slate. So I started by completely re-installing mambaforge. Since
my tmux is in mambaforge, I had to first make sure all my tmux sessions
were killed.

``` sh
(base) [login10: ~]--% which tmux
/projects/eriq@colostate.edu/mambaforge/bin/tmux
(base) [login10: ~]--% tmux kill-server
error connecting to /tmp/tmux-2002325/default (No such file or directory)

# It looks like nothing was actually running
# Wait! That is because I typically login to login11...
# as per my .ssh/config:

Host summit
   HostName login11.rc.colorado.edu
   RemoteForward 52611 localhost:52698

# So, I got onto that instead.
ssh eriq@colostate.edu@summit

# but I had closed all my tmux windows already.  Good
(base) [login11: ~]--% tmux kill-server
no server running on /tmp/tmux-2002325/default
```

Now, it is time to reinstall mambaforge.

First, I removed the conda block from my .bashrc, logged out and logged
back in again.

Then I removed my mambaforge directory:

``` sh
[login11: projects]--% pwd
/home/eriq@colostate.edu/projects
[login11: projects]--% ls
afblue  Cons-Gen-Comp-2022  DaddySalmon  empty_dir  java-programs  mambaforge  README.mdwn  rmote-server
[login11: projects]--% rm -rf mambaforge
```

That takes a little bit of time. While it was running, I was reading
through the Alpine docs.

Once it was done, I went ahead and did a fresh install of mambaforge as
detailed here:
<https://mamba.readthedocs.io/en/latest/installation.html>

``` sh
[login11: projects]--% acompile
acompile: submitting job... salloc --nodes=1 --partition=acompile --ntasks=1  --time=01:00:00  --qos=compile --job-name=acompile --bell --oversubscribe srun --pty /bin/bash
salloc: Granted job allocation 2291926
salloc: Nodes c3cpu-c11-u21-2 are ready for job

wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

# after agreeing to the licensing, I had to tell it to put the thing in:
/projects/eriq@colostate.edu/mambaforge

# while it was installing, I noticed that it is using python 3.10
Extracting python-3.10.12-hd12c33a_0_cpython.conda
```

That was painless and incredibly fast. After it was done, I sourced my
.bashrc to initialize it:

``` sh
[c3cpu-c11-u21-2: projects]--% source ~/.bashrc
(base) [c3cpu-c11-u21-2: projects]--%
```

Once in my base environment, I installed tmux there. I simply installed
the latest tmux on conda: 3.5

``` sh
(base) [c3cpu-c11-u21-2: projects]--% mamba install -c conda-forge tmux
```

Once that was done. I logged out and logged back in via my iTerm tmux.

# Get Snakemake

I have recently installed and tested snakemake 7.30.1 on SEDNA and it is
working great. That is what is still currently up on Anaconda so I will
use it here, as well, for consistency. Following the instructions at:
<https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>
I did:

``` sh
(base) [login11: ~]--% acompile
acompile: submitting job... salloc --nodes=1 --partition=acompile --ntasks=1  --time=01:00:00  --qos=compile --job-name=acompile --bell --oversubscribe srun --pty /bin/bash
salloc: Granted job allocation 2291932
salloc: Nodes c3cpu-c11-u21-2 are ready for job
(base) [c3cpu-c11-u21-2: ~]--% mamba create -c conda-forge -c bioconda -n snakemake-7.30.1  snakemake
```

That was painless.

# Get mega-non-model-wgs snakeflow and test it in an interactive session

I symlinked `scratch` in my home directory to alpine’s scratch and am
now there.

``` sh
(base) [login11: ~]--% cd scratch/
(base) [login11: scratch]--% ls
README.mdwn
(base) [login11: scratch]--% gitup
Agent pid 31431
Enter passphrase for /home/eriq@colostate.edu/.ssh/id_ed25519:
Identity added: /home/eriq@colostate.edu/.ssh/id_ed25519 (SUMMIT)
(base) [login11: scratch]--%
(base) [login11: scratch]--% git clone git@github.com:eriqande/mega-non-model-wgs-snakeflow.git
```

The `gitup` there is an alias I have to refresh my SSH credentials,
which seems necessary occasionally on Alpine. The alias looks like this
in my .bashrc:

``` sh
# here is an alias for the command you need to do to
# make SUMMIT have access to your SSH keys for GitHub.
# It seems like I have to do this pretty much every new
# session, and sometime after a certain amount of time
# in the same session, so I just made it an alias. Also,
# it requires my SSH passphrase, so I have to do it
# interactively.
alias gitup='eval "$(ssh-agent -s)"; ssh-add ~/.ssh/id_ed25519'
```

At any rate, now I will first get an interactive session with a single
core to first install all of the conda environments for the workflow:

``` sh
(base) [login11: scratch]--% cd mega-non-model-wgs-snakeflow/
(base) [login11: mega-non-model-wgs-snakeflow]--% module load slurm/alpine
(base) [login11: mega-non-model-wgs-snakeflow]--% sinteractive -t 00:20:00 -c 1
/usr/local/bin/sinteractive: waiting for job with ID 2292145 to start.
```

I note that using `sinteractive` gets `screen` involved, and I am not a
big fan of that when I am already using tmux. So, let’s see if srun will
work for us:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% srun -c 1 -t 00:20:00 --pty /bin/bash
```

Here on the fourth of july, 2023 in the early morning, it did not take
long to get those resources!

``` sh
(base) [c3cpu-c9-u1-1: mega-non-model-wgs-snakeflow]--% conda activate snakemake-7.30.1
(snakemake-7.30.1) [c3cpu-c9-u1-1: mega-non-model-wgs-snakeflow]--% snakemake --cores 1 --use-conda --conda-create-envs-only --configfile .test/config/config.yaml
```

We will see if all of those can get installed in less than 20 minutes
before my allocation dies!

That seemed to work. So, now I will get 4 cores on atesting and see how
it works:

``` sh
srun --partition=atesting -c 4 -t 00:30:00 --pty /bin/bash

# dry-run
snakemake --cores 4 -np --use-conda --configfile .test/config/config.yaml

# that looked fine, so I did a full fun
```

That ran without any problems on the little test case. 322 steps done
with no problems at all.

So, what we have left to do now is try it using the slurm profile on a
real data set, like the SWTHs.

# Start copying data over

My token needed refreshing. I ended up just making a new rclone remote
for my onedrive and then put a link to the RueggLab_BGPshare into my
onedrive root directory.

I am going to copy using two different commands in different shells.
First is:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% fp .
/scratch/alpine/eriq@colostate.edu/mega-non-model-wgs-snakeflow/.
(base) [login11: mega-non-model-wgs-snakeflow]--% rclone copy onedrive2:BGP_Share/Genetic_and_Environmental_Data/Novoseq_fastq_2_download/SWTH  data/SWTH

# Holy freakazoid that is fast!  About 400 Mb/sec
```

The second copy job is:

``` sh
rclone copy  onedrive2:BGP_Share/Genetic_and_Environmental_Data/Novoseq_fastq_2_download/SWTH_Novoseq_Plate2 data/SWTH_Novoseq_Plate2
```

# copy the old resources directory that has the indexed genome

This is up on the sharepoint:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% pwd
/home/eriq@colostate.edu/scratch/mega-non-model-wgs-snakeflow
(base) [login11: mega-non-model-wgs-snakeflow]--% rm -r resources

(base) [login11: mega-non-model-wgs-snakeflow]--% rclone copy  onedrive2:BGP_Share/Genetic_and_Environmental_Data/Species_genetic_data/SWTH/resources  resources
```

# Make an alpine slurm profile

We have to change the partitions, etc. But more importantly we will need
to modify the memory requirements and the tmpdir in resources and maybe
the threads.

We are going to do that in the profile.

We do some exploring first. Here are all the places where the resources
get modified in the various rules:

``` sh
awk 'BEGIN {OFS="\t"} /^rule/ {rule=$2; next} /resources:/ {go=1; print FILENAME, rule, $0;  next} /[a-z]+:/ {go=0} /=/ && go==1 {print FILENAME, rule, $0}' ../Snakefile *.smk
calling.smk make_gvcf_sections:     resources:
calling.smk make_gvcf_sections:         time="1-00:00:00",
calling.smk make_gvcf_sections:         mem_mb = 4600,
calling.smk make_gvcf_sections:         cpus = 1
calling.smk genomics_db_import_chromosomes:     resources:
calling.smk genomics_db_import_chromosomes:         mem_mb = 9400,
calling.smk genomics_db_import_chromosomes:         cpus = 2,
calling.smk genomics_db_import_chromosomes:         time = "36:00:00"
calling.smk genomics_db_import_chromosomes: # job.  resources: cpus = 2, does not work.  Same for above. I set reader
calling.smk genomics_db_import_scaffold_groups:     resources:
calling.smk genomics_db_import_scaffold_groups:         mem_mb = 9400,
calling.smk genomics_db_import_scaffold_groups:         cpus = 2,
calling.smk genomics_db_import_scaffold_groups:         time = "36:00:00"
calling.smk genomics_db2vcf_scattered:      resources:
calling.smk genomics_db2vcf_scattered:          mem_mb = 11750,
calling.smk genomics_db2vcf_scattered:          cpus = 2,
calling.smk genomics_db2vcf_scattered:          time = "1-00:00:00"
mapping.smk mark_duplicates:        resources:
mapping.smk mark_duplicates:            cpus = 1
qc.smk  multiqc_dir:        resources:
qc.smk  multiqc_dir:            mem_mb = 36800
ref.smk bwa_index:      resources:
ref.smk bwa_index:          mem_mb=36900,
```

And here are the only places where we are using more than 1 thread:

``` sh
awk 'BEGIN {OFS="\t"} /^rule/ {rule=$2; next} /threads:/ {go=1; print FILENAME, rule, $0;  next}' ../Snakefile *.smk
angsd-ready-bams.smk    realigner_target_creator:       threads: 4
annotation.smk  annotate_variants:      threads: 4
calling.smk make_gvcf_sections:     threads: 1
calling.smk genomics_db_import_chromosomes:     threads: 2
calling.smk genomics_db_import_scaffold_groups:     threads: 2
calling.smk genomics_db2vcf_scattered:      threads: 2
mapping.smk map_reads:      threads: 4
```

All the nodes we would send jobs to on Alpine have 3.74 Gb of RAM per
core. So we merely need to make sure that the threads and the memory are
in agreement and then off we go…

I modified the profile accordingly.

# Running

I ran this first on SWTH from a login node with:

``` sh
snakemake --profile hpcc-profiles/slurm/alpine --configfile config/config.yaml
```

It worked better than expected. But:

- I had allowed too many jobs/cores and it submitted more than was
  allowed. So some jobs were killed right after submission.
- I also got a lot of OOM error on the map_reads rule. The problem there
  is that the resources don’t seem to be getting set properly by the
  –profile values.

I will be able to explore that some more. I will try putting them into
the –workflow-profile

To kill it gracefully, from another terminal I grepped snakemake out of
ps -ax and found the job id number and then just used `kill` to
gracefully shut it down. It continued to run the jobs that were running
(and that were queued…but I scancelled those) until those things
finished. All in all, a very nice experience.

There is not good documentation on the max number of jobs, but I will
scale it down to 600 in the profile.

Here are some of the samples that failed on map_reads:

``` sh
(base) [login11: map_reads]--% grep oom * | awk -F"," '{print $2}' | head -n 20
sample=00N0540
sample=00N2156
sample=00N2157
sample=00N2158
sample=00N2159
sample=00N2160
sample=00N2161
sample=00N5917
sample=00N5929
sample=07N31126
sample=07N31128
sample=07N51683
sample=07N51685
sample=07N51690
sample=07N51755
sample=07N51761
sample=07N51764
sample=07N51768
sample=07N51770
sample=07N51771
```

So, I can test a few of these when I crank the memory back up.

*It might be that it needs fewer cores but more memory* The mem might
scale with the number of threads…

I figured out what was wrong with the profile and it is now working much
better.

Also, here is how to compute the total number of CPU hours (SU’s) used
from a certain date:

``` sh
sacct -S2023-01-01 -u eriq@colostate.edu -ojobid,start,end,alloccpu,cputime | column -t | awk '!/[a-z.]/ {n=split($5,a,/:/); h=a[1]; m=a[2]; s=a[3]; fh = h + m/60 + s/3600; cumul += $4 * fh; print $0, fh} END {print "TotalCoreHours:", cumul}'
```

Right now I am here:

``` sh
2318029       2023-07-06T15:08:31  2023-07-06T15:13:14  4           00:18:52 0.314444
2318033       2023-07-06T15:09:46  2023-07-06T15:15:25  4           00:22:36 0.376667
2318042       2023-07-06T15:09:46  2023-07-06T15:14:49  4           00:20:12 0.336667
2318059       2023-07-06T15:09:46  2023-07-06T15:13:35  4           00:15:16 0.254444
2318069       2023-07-06T15:10:30  2023-07-06T15:14:32  4           00:16:08 0.268889
2318095       2023-07-06T15:10:30  2023-07-06T15:12:42  4           00:08:48 0.146667
2318126       2023-07-06T15:11:47  2023-07-06T15:15:57  4           00:16:40 0.277778
2318133       2023-07-06T15:12:40  2023-07-06T15:15:26  4           00:11:04 0.184444
TotalCoreHours: 312.279
```

OK!

# Restarting

That thing ran for about 12 hours and got quite far, but then ended.
There had been some errors.

The dry-run says there are 5 mark_duplicates and 7 clip_overlaps to be
done, and then all the other stuff downstream of that.

By grepping out of the log file that snakemake automatically saves we
can get the slurm job numbers that failed:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% grep Error .snakemake/log/2023-07-06T145637.260179.snakemake.log | awk '/executing/ {err[$4] = sprintf("%s %s", err[$4], $10)} END {for(i in err) print i, err[i]}'
make_gvcf_sections  2325409 2325410 2325411 2325412 2325413 2325414 2325415 2325416 2325417 2325419 2325420 2325421 2325422 2325423 2325424 2325425 2325426 2325427 2325428 2325429 2325430
clip_overlaps  2325408 2325418
mark_duplicates  2321682 2322604 2330469 2337276 2337861
```

And now we can investigate those failures at our leisure.

Here is another way to print those that makes it easier to use them:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% grep Error .snakemake/log/2023-07-06T145637.260179.snakemake.log | awk '/executing/ {err[$4] = sprintf("%s *%s*", err[$4], $10)} END {for(i in err) print i, err[i]}' | sed 's/,//g;'
make_gvcf_sections  *2325409* *2325410* *2325411* *2325412* *2325413* *2325414* *2325415* *2325416* *2325417* *2325419* *2325420* *2325421* *2325422* *2325423* *2325424* *2325425* *2325426* *2325427* *2325428* *2325429* *2325430*
clip_overlaps  *2325408* *2325418*
mark_duplicates  *2321682* *2322604* *2330469* *2337276* *2337861*
```

Note that some of the make_gvcf_sections failed.

``` sh
# explore these:
# The mark_duplicates fails have slurm logs:
(base) [login11: mark_duplicates]--% ls -l *2321682* *2322604* *2330469* *2337276* *2337861*
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 4021 Jul  6 17:27 mark_duplicates-bqsr_round=0,sample=00N5507-2322604.err
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu    0 Jul  6 17:26 mark_duplicates-bqsr_round=0,sample=00N5507-2322604.out
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 3986 Jul  6 21:08 mark_duplicates-bqsr_round=0,sample=08N0723-2337861.err
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu    0 Jul  6 21:02 mark_duplicates-bqsr_round=0,sample=08N0723-2337861.out
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 4045 Jul  6 17:05 mark_duplicates-bqsr_round=0,sample=13N01565-2321682.err
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu    0 Jul  6 17:03 mark_duplicates-bqsr_round=0,sample=13N01565-2321682.out
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 3986 Jul  6 20:28 mark_duplicates-bqsr_round=0,sample=99N0667-2337276.err
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu    0 Jul  6 20:24 mark_duplicates-bqsr_round=0,sample=99N0667-2337276.out
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 3986 Jul  6 18:55 mark_duplicates-bqsr_round=0,sample=99N0674-2330469.err
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu    0 Jul  6 18:51 mark_duplicates-bqsr_round=0,sample=99N0674-2330469.out

# and, in every case they were out-of-memory kills:
(base) [login11: mark_duplicates]--% grep kill  *2321682* *2322604* *2330469* *2337276* *2337861*
mark_duplicates-bqsr_round=0,sample=13N01565-2321682.err:slurmstepd: error: Detected 1 oom-kill event(s) in StepId=2321682.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.
mark_duplicates-bqsr_round=0,sample=00N5507-2322604.err:slurmstepd: error: Detected 1 oom-kill event(s) in StepId=2322604.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.
mark_duplicates-bqsr_round=0,sample=99N0674-2330469.err:slurmstepd: error: Detected 1 oom_kill event in StepId=2330469.batch. Some of the step tasks have been OOM Killed.
mark_duplicates-bqsr_round=0,sample=99N0667-2337276.err:slurmstepd: error: Detected 1 oom_kill event in StepId=2337276.batch. Some of the step tasks have been OOM Killed.
mark_duplicates-bqsr_round=0,sample=08N0723-2337861.err:slurmstepd: error: Detected 1 oom_kill event in StepId=2337861.batch. Some of the step tasks have been OOM Killed.
```

But I don’t have slurm logs for the other failures. So they must have
failed before they got launched by SLURM. Weird…

Let’s see if we can find the logs for them:

``` sh
# here are the logs for them:
(base) [login11: mega-non-model-wgs-snakeflow]--% grep -A 10  'Error in rule make_gvcf_sections'  .snakemake/log/2023-07-06T145637.260179.snakemake.log | grep log:
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046221.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046221.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046224.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046224.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/07N51764/NC_046223.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/07N51764/NC_046223.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/00N5500/NC_046258.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/00N5500/NC_046258.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046221.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046221.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046225.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046225.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/07N52116/NC_046233.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/07N52116/NC_046233.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/07N51764/NC_046226.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/07N51764/NC_046226.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046233.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046233.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046224.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046224.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046233.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046233.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/05N5310/NC_046252.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/05N5310/NC_046252.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046235.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/AB16098/NC_046235.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/14N02035/NC_046230.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/14N02035/NC_046230.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046230.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/97N6564/NC_046230.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/96N0608/NC_046247.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/96N0608/NC_046247.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/98N1131/NC_046262.2.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/98N1131/NC_046262.2.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/99N0668/NC_046251.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/99N0668/NC_046251.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/00N5502/NC_046225.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/00N5502/NC_046225.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/07N53531/NC_046230.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/07N53531/NC_046230.1.stdout (check log file(s) for error details)
    log: results/bqsr-round-0/logs/gatk/haplotypecaller/21N01277/NC_046246.1.stderr, results/bqsr-round-0/logs/gatk/haplotypecaller/21N01277/NC_046246.1.stdout (check log file(s) for error details)
```

And those are all empty.

So, now let us check sacct for those job IDs.

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2325409,2325410,2325411
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
2325409      smk-make_+     amilan csu-gener+          1     FAILED     0:53
2325409.bat+      batch            csu-gener+          1  CANCELLED     0:53
2325409.ext+     extern            csu-gener+          1  COMPLETED      0:0
2325410      smk-make_+     amilan csu-gener+          1     FAILED     0:53
2325410.bat+      batch            csu-gener+          1  CANCELLED     0:53
2325410.ext+     extern            csu-gener+          1  COMPLETED      0:0
2325411      smk-make_+     amilan csu-gener+          1     FAILED     0:53
2325411.bat+      batch            csu-gener+          1  CANCELLED     0:53
2325411.ext+     extern            csu-gener+          1  COMPLETED      0:0
```

Hey! Those jobs just got cancelled before running. Probably no big deal.

And the clip-overlaps:

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2325408,2325418
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
2325408      smk-clip_+     amilan csu-gener+          1     FAILED     0:53
2325408.bat+      batch            csu-gener+          1  CANCELLED     0:53
2325408.ext+     extern            csu-gener+          1  COMPLETED      0:0
2325418      smk-clip_+     amilan csu-gener+          1     FAILED     0:53
2325418.bat+      batch            csu-gener+          1  CANCELLED     0:53
2325418.ext+     extern            csu-gener+          1  COMPLETED      0:0
```

Same thing—just got cancelled.

And finally, let’s check the other mark_duplicate runs this way. Here
are the numbers in an easier to use format:

``` sh
grep Error .snakemake/log/2023-07-06T145637.260179.snakemake.log | awk '/executing/ {err[$4] = sprintf("%s%s", err[$4], $10)} END {for(i in err) print i, err[i]}'
make_gvcf_sections 2325409,2325410,2325411,2325412,2325413,2325414,2325415,2325416,2325417,2325419,2325420,2325421,2325422,2325423,2325424,2325425,2325426,2325427,2325428,2325429,2325430,
clip_overlaps 2325408,2325418,
mark_duplicates 2321682,2322604,2330469,2337276,2337861
```

``` sh
(base) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2321682,2322604,2330469,2337276,2337861
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
2321682      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2321682.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2321682.ext+     extern            csu-gener+          1 OUT_OF_ME+    0:125
2322604      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2322604.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2322604.ext+     extern            csu-gener+          1 OUT_OF_ME+    0:125
2330469      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2330469.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2330469.ext+     extern            csu-gener+          1  COMPLETED      0:0
2337276      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2337276.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2337276.ext+     extern            csu-gener+          1  COMPLETED      0:0
2337861      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2337861.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2337861.ext+     extern            csu-gener+          1  COMPLETED      0:0
```

Yep! All out of memories. So we just give those more memory when
starting it all up again. It does not look like it ran out of Java
memory—the SLURM controller killed it. So we will just restart it with
three times as much memory.

It only has the default 3700 Mb so, let’s give it 11220 Mb. We can do
that on the command line with –set-resources:

``` sh
# check it with a dry run
snakemake  -np  --profile hpcc-profiles/slurm/alpine --configfile config/config.yaml --set-resources mark_duplicates:mem_mb=11220

# that dry run showed that the memory had been pumped up on the mark duplicates runs
# so, now, start it up for real.
snakemake   --profile hpcc-profiles/slurm/alpine --configfile config/config.yaml --set-resources mark_duplicates:mem_mb=11220
```

That run finished after a while, but three of the mark_duplicates jobs
ran out of memory (and these are the biggest files, too). So, I am going
to start them again and give them 8 cores worth of memory: 29920.

``` sh
snakemake   --profile hpcc-profiles/slurm/alpine --configfile config/config.yaml --set-resources mark_duplicates:mem_mb=29920
```

Hey! Check this out…if you apply `--set-resources` on the command line,
it completely overwrites the values set with `set-resources:` in the
profile. That is kind of BS, I’ve gotta say.

But, good to know. The take-home message is: if you want to change the
value of one resource for one rule, just go ahead and make that change
in the profile, not on the command line. Lame! Sad! Bogus! (But I guess
I can see why that is…multiple invocations of –set-resources completely
overwrite previous ones…) As a consequence of this choice, however, I am
going to have to restart this because the import_genomics_db rules
reverted back to the 36 hour time limit, which exceeds the QoS. No
worries. I soft-killed it and am waiting for the final make_gvcf rules
to finish out.

When it was done I just started it up again with:

``` sh
snakemake   --profile hpcc-profiles/slurm/alpine --configfile config/config.yaml
```

That ran along fine. Until I got an error:

``` sh
sacct failed to return the status for jobid 2355542
Maybe you need to use scontrol instead?
Failed to obtain job status. See above for error message.
```

Checking it myself I see this:

``` sh
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2355542 --parsable2
JobID|JobName|Partition|Account|AllocCPUS|State|ExitCode
2355542|smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0020,sg_or_chrom=NC_046222.1|amilan|csu-general|2|FAILED|1:0
2355542.batch|batch||csu-general|2|FAILED|1:0
2355542.extern|extern||csu-general|2|COMPLETED|0:0

# and here is what the slurm log says:
(base) [login11: genomics_db2vcf_scattered]--% head -n 200  *2355542*
==> genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0020,sg_or_chrom=NC_046222.1-2355542.err <==
Waiting at most 60 seconds for missing files.
Traceback (most recent call last):
  File "<frozen runpy>", line 198, in _run_module_as_main
  File "<frozen runpy>", line 88, in _run_code
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/site-packages/snakemake/__main__.py", line 4, in <module>
    main()
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/site-packages/snakemake/__init__.py", line 3079, in main
    wait_for_files([args.wait_for_files_file], latency_wait=args.latency_wait)
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/site-packages/snakemake/io.py", line 878, in wait_for_files
    raise IOError(
OSError: Missing files after 60 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
/scratch/alpine/eriq@colostate.edu/mega-non-model-wgs-snakeflow/.snakemake/tmp.wxmfac_p/snakejob.genomics_db2vcf_scattered.18166.sh.waitforfilesfile.txt
```

The main job does not exist. Whoa! This job did not create a log at all!
It is like it died super early on. Perhaps the genomics db is corrupted.

In fact all of the genomics_db2vcf_scattered rule runs for
sg_or_chrom=NC_046222.1 failed in exactly the same way. Perhaps the
genomicsdbimport rule was done but the files had not all appeared,
because of filesystem latency, but the receipts were there, and so all
the genotypeGVCFs jobs got confused and failed as soon as they were
launched (they were all launched with 2 seconds of one another).

These jobs are still running, even though snakemake aborted itself:

``` sh
       JOBID PARTITION                                                                                                 NAME       USER ST            TIME  NODES   NODELIST(REASON)  CPUS   MIN_MEMORY   TIME_LIMIT    TIME_LEFT PRIORITY
     2355521    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0008,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355522    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0037,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355523    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0030,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355524    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0023,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355525    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0051,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355526    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0020,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355527    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0048,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355528    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0041,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355529    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0006,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355530    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0034,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355531    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0027,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355532    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0055,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355533    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0045,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355534    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0017,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355535    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0009,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355536    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0003,sg_or_chrom=NC_046221.1 eriq@colos  R           25:40      1    c3cpu-c11-u34-4     2       11000M   1-00:00:00     23:34:20 0.00000378931873
     2355501    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0043,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a9-u5-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355502    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0038,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a9-u5-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355503    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0013,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a9-u5-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355504    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0036,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1     c3cpu-a5-u15-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355505    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0001,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1     c3cpu-a5-u15-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355506    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0029,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1     c3cpu-a5-u15-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355507    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0022,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1     c3cpu-a5-u15-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355508    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0050,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1     c3cpu-a5-u15-2     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355509    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0019,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355510    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0047,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355511    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0011,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355512    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0040,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355513    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0007,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355514    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0005,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355515    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0033,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355516    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0002,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355517    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0026,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355518    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0054,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355519    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0016,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355520    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0044,sg_or_chrom=NC_046221.1 eriq@colos  R           26:10      1      c3cpu-a2-u3-1     2       11000M   1-00:00:00     23:33:50 0.00000378931873
     2355481    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0024,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-c9-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355482    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0052,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-c9-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355483    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0031,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1     c3cpu-c11-u3-3     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355484    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0014,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1     c3cpu-c11-u3-3     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355485    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0042,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1     c3cpu-c11-u3-4     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355486    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0035,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a9-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355487    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0028,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a9-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355488    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0056,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a7-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355489    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0021,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a7-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355490    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0049,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a7-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355491    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0018,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a7-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355492    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0046,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a5-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355493    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0039,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a5-u1-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355494    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0010,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a2-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355495    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0004,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a2-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355496    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0032,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a2-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355497    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0012,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-a2-u1-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355498    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0025,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1      c3cpu-c9-u5-2     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355499    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0053,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1     c3cpu-c11-u5-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355500    amilan                 smk-genomics_db2vcf_scattered-bqsr_round=0,scatter=scat_0015,sg_or_chrom=NC_046221.1 eriq@colos  R           26:40      1     c3cpu-c11-u5-1     2       11000M   1-00:00:00     23:33:20 0.00000378931873
     2355540    amilan                            smk-bung_filtered_vcfs_back_together-bqsr_round=0,sg_or_chrom=NC_046224.1 eriq@colos  R           23:09      1      c3cpu-c9-u1-1     1        3740M      8:00:00      7:36:51 0.00000378908590
     2355472    amilan                                         smk-mark_dp0_as_missing-bqsr_round=0,sg_or_chrom=NC_046223.1 eriq@colos  R           35:43      1      c3cpu-c9-u1-1     1        3740M      8:00:00      7:24:17 0.00000378908590
```

I will check the logs on these when they are done to make sure that they
finished correctly.

I did check that. They were all fine. I restarted it and then we ran
into two out-of-memories on:

    (base) [login11: mega-non-model-wgs-snakeflow]--% grep  Error  .snakemake/log/2023-07-08T173040.213753.snakemake.log
    Error in rule genomics_db2vcf_scattered:
    Error executing rule genomics_db2vcf_scattered on cluster (jobid: 19130, external: 2357623, jobscript: /scratch/alpine/eriq@colostate.edu/mega-non-model-wgs-snakeflow/.snakemake/tmp.ybs08jv_/snakejob.genomics_db2vcf_scattered.19130.sh). For error details see the cluster log and the log files of the involved rule(s).
    Error in rule genomics_db2vcf_scattered:
    Error executing rule genomics_db2vcf_scattered on cluster (jobid: 19128, external: 2357643, jobscript: /scratch/alpine/eriq@colostate.edu/mega-non-model-wgs-snakeflow/.snakemake/tmp.ybs08jv_/snakejob.genomics_db2vcf_scattered.19128.sh). For error details see the cluster log and the log files of the involved rule(s).

Predictably, these were scat006 and scat007 on scaff_group_001. We can
check the logs to see how far they got.

They only got about half way through. Such BS. I am just going to give
each sequence its own scatter group in those and let it work it out that
way:

``` sh
(base) [login11: config]--% pwd
/home/eriq@colostate.edu/scratch/mega-non-model-wgs-snakeflow/config

(base) [login11: config]--% awk 'BEGIN {OFS="\t"} $1=="scaff_group_001" && $2~/scat_000[67]/ {new2 = $2 "." ++n; print $1, new2, $3, $4, $5, $6;  next} {print}' scatters_3000000.tsv > scatters_3000000_exploded.tsv

# and then we edit the config.yaml to use scatters_3000000_exploded.tsv
```

After doing a dry run, it was clear that we need to set the modification
date for scatters_3000000_exploded.tsv to earlier since that file is
considered an input, and so it will want to rerun all the GenotypeGvcfs
steps. No worries:

``` sh
touch -d "last Monday" config/scatters_3000000_exploded.tsv

# and check:
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% ls -l config/scatters_3000000_exploded.tsv
-rw-r--r--. 1 eriq@colostate.edu eriqgrp@colostate.edu 30371 Jul  3 00:00 config/scatters_3000000_exploded.tsv
```

(Pretty rad that you can use “last Monday” to do that!)

And that finished up by Sunday noon. So 3-days total to get that
processed. Not bad.

And check how many CPU hours that was:

``` sh
sacct -S2023-07-03 -u eriq@colostate.edu -ojobid,start,end,alloccpu,cputime | column -t | awk '!/[a-z.]/ {n=split($5,a,/:/); h=a[1]; m=a[2]; s=a[3]; fh = h + m/60 + s/3600; cumul += $4 * fh;} END {print "TotalCoreHours:", cumul}'
TotalCoreHours: 8875.77
```

OK!

# Reboot! The sample IDs were not the same between runs

There were about 70 birds that had been re-run, but the seq center put a
Z in front of their names, so I didn’t realize they were the same birds.
So, the units file will need to be re-made to show two units for each of
those birds.

``` r
# doing this on my laptop
library(tidyverse)

old_units <- read_tsv("~/Documents/git-repos/mega-non-model-configs/SWTH-2023-07-05/config/OLD_units.tsv")

# learn a bit about these.  There are 74 samples starting with a Z
old_units %>%
  count(str_detect(sample, "^Z"))
```

So, which of those are re-runs?

``` r
old_units %>%
  mutate(noZ = str_remove(sample, "^Z")) %>%
  count(noZ) %>%
  count(n)
```

So, all 74 of those are re-runs. Thus, we can just take the Z off the
sample and sample_id fields and then put them together and amend the
units column:

``` r
new_units <- old_units %>%
  mutate(
    sample = str_remove(sample, "^Z"),
    sample_id = str_remove(sample_id, "^Z")
  ) %>%
  arrange(sample, library) %>%
  group_by(sample) %>%
  mutate(unit = 1:n())
```

But, I need to ask CH if these were completely new library preps or not,
because that affects the duplicate marking stage. OK! They were
completely new library preps, so new_units is correct as it it, with
different library preps in there for the resequences.

So, now I can save that as the units file

``` r
write_tsv(new_units, file = "~/Documents/git-repos/mega-non-model-configs/SWTH-2023-07-05/config/units.tsv")
```

Then, when I used that new units file with:

``` sh
snakemake -np --profile  hpcc-profiles/slurm/alpine  --configfile config/config.yaml
```

it said everything was done. I am a little surprised by that, because it
seems like it should notice that the inputs are different to some of the
rules. But perhaps not, since the file modification times are not any
different on all those. So, there are various ways that I could conceive
of dealing with that…

For fun, I could try touching the data files for the 74 “Z” fastqs. But
the problem with that is that it will trigger new fastqc runs that are
not necessary. But, that might be the cleanest way to do it anyway. I
will try:

``` sh
(base) [login11: SWTH_Novoseq_Plate2]--% pwd
/home/eriq@colostate.edu/scratch/mega-non-model-wgs-snakeflow/data/SWTH_Novoseq_Plate2
(base) [login11: SWTH_Novoseq_Plate2]--% touch  Z*/*fq.gz

# that was two files for each of the 74 reruns:
(base) [login11: SWTH_Novoseq_Plate2]--% ls -l   Z*/*fq.gz | wc
    148    1332   20164
```

Now, when I give it a dry run, I see:

``` sh
Job stats:
job                     count    min threads    max threads
--------------------  -------  -------------  -------------
all                         1              1              1
clip_overlaps              74              1              1
concat_gvcf_sections       74              1              1
fastqc_read1               74              1              1
fastqc_read2               74              1              1
make_gvcf_sections       3108              1              1
map_reads                 148              4              4
mark_duplicates            74              1              1
multiqc_dir                 1              1              1
samtools_stats             74              1              1
trim_reads_pe             148              1              1
total                    3850              1              4
```

Note that there are 41 “chromosomes” and 1 scaffold group. And 74 \* 42
= 3108.

So, that looks to be right. It is going to have to re-trim everything,
because I deleted the trimmed files, but that is not such a big deal.

Also, those fastqc’s do have to be redone because the naming of them
have changed.

But, what is interesting here is that this has not triggered a re-run of
the genomics data bases and beyond. I expect that is because of the way
that I rigged things so that I could add individuals to the genomics
data bases. At some point I should jettison that stuff. But for now I
can just delete all the genomicsDBs. They have to be redone anyway!

``` sh
(base) [login11: bqsr-round-0]--% pwd
/home/eriq@colostate.edu/scratch/mega-non-model-wgs-snakeflow/results/bqsr-round-0
(base) [login11: bqsr-round-0]--% rm -rf gdb_intervals gdb_accounting genomics_db
```

Now that I have done that, my dry run looks like:

``` sh
Job stats:
job                                   count    min threads    max threads
----------------------------------  -------  -------------  -------------
all                                       1              1              1
bcf_concat                                1              1              1
bcf_concat_mafs                           1              1              1
bcf_maf_section_summaries                42              1              1
bcf_section_summaries                   126              1              1
bung_filtered_vcfs_back_together         42              1              1
clip_overlaps                            74              1              1
combine_bcftools_stats                    3              1              1
combine_maf_bcftools_stats                1              1              1
concat_gvcf_sections                     74              1              1
fastqc_read1                             74              1              1
fastqc_read2                             74              1              1
gather_scattered_vcfs                    42              1              1
genomics_db2vcf_scattered               494              2              2
genomics_db_import_chromosomes           41              2              2
genomics_db_import_scaffold_groups        1              2              2
hard_filter_indels                       42              1              1
hard_filter_snps                         42              1              1
maf_filter                               42              1              1
make_gvcf_sections                     3108              1              1
make_indel_vcf                           42              1              1
make_snp_vcf                             42              1              1
map_reads                               148              4              4
mark_dp0_as_missing                      42              1              1
mark_duplicates                          74              1              1
multiqc_dir                               1              1              1
samtools_stats                           74              1              1
trim_reads_pe                           148              1              1
total                                  4896              1              4
```

That looks just about right! It also looks like I have removed all of
the `protected()` statements in the workflow, so this should not run
into any problems. Let’s load this bad boy up!

``` sh
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% snakemake -p --profile  hpcc-profiles/slurm/alpine  --configfile config/config.yaml
```

By morning that had ended with some failures. It was easy to find them:

``` sh
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% grep Error .snakemake/log/2023-07-27T160926.747195.snakemake.log  | awk '/executing/ {err[$4] = sprintf("%s%s", err[$4], $10)} END {for(i in err) print i, err[i]}'

mark_duplicates 2516965,2517268,2517263,2517378,2517370,2518710
```

And those were out of memory, and some that just failed.

``` sh
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2516965,2517268,2517263,2517378,2517370,2518710
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
2516965      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2516965.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2516965.ext+     extern            csu-gener+          1  COMPLETED      0:0
2517263      smk-mark_+     amilan csu-gener+          1     FAILED      1:0
2517263.bat+      batch            csu-gener+          1     FAILED      1:0
2517263.ext+     extern            csu-gener+          1  COMPLETED      0:0
2517268      smk-mark_+     amilan csu-gener+          1     FAILED      1:0
2517268.bat+      batch            csu-gener+          1     FAILED      1:0
2517268.ext+     extern            csu-gener+          1  COMPLETED      0:0
2517370      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2517370.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2517370.ext+     extern            csu-gener+          1  COMPLETED      0:0
2517378      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2517378.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2517378.ext+     extern            csu-gener+          1  COMPLETED      0:0
2518710      smk-mark_+     amilan csu-gener+          1 OUT_OF_ME+    0:125
2518710.bat+      batch            csu-gener+          1 OUT_OF_ME+    0:125
2518710.ext+     extern            csu-gener+          1  COMPLETED      0:0
```

So, I will just change the profile to give mark_duplicates more memory
and then I will restart.

Here is what remains:

``` sh
Job stats:
job                                   count    min threads    max threads
----------------------------------  -------  -------------  -------------
all                                       1              1              1
bcf_concat                                1              1              1
bcf_concat_mafs                           1              1              1
bcf_maf_section_summaries                42              1              1
bcf_section_summaries                   126              1              1
bung_filtered_vcfs_back_together         42              1              1
clip_overlaps                             6              1              1
combine_bcftools_stats                    3              1              1
combine_maf_bcftools_stats                1              1              1
concat_gvcf_sections                      6              1              1
gather_scattered_vcfs                    42              1              1
genomics_db2vcf_scattered               494              2              2
genomics_db_import_chromosomes           41              2              2
genomics_db_import_scaffold_groups        1              2              2
hard_filter_indels                       42              1              1
hard_filter_snps                         42              1              1
maf_filter                               42              1              1
make_gvcf_sections                      252              1              1
make_indel_vcf                           42              1              1
make_snp_vcf                             42              1              1
mark_dp0_as_missing                      42              1              1
mark_duplicates                           6              1              1
multiqc_dir                               1              1              1
samtools_stats                            6              1              1
total                                  1324              1              2
```

It seems to be working great.

It did as much as it could but it threw two errors on genomics_db2vcf
rule runs. sacct says that both of the jobs finished, so it must just be
a job latency thing.

It said:

``` sh
sacct failed to return the status for jobid 2527084
Maybe you need to use scontrol instead?
Failed to obtain job status. See above for error message.
WorkflowError:
Failed to obtain job status. See above for error message.
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/asyncio/runners.py", line 190, in run
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/asyncio/runners.py", line 118, in run
  File "/projects/eriq@colostate.edu/mambaforge/envs/snakemake-7.30.1/lib/python3.11/asyncio/base_events.py", line 653, in run_until_complete
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: .snakemake/log/2023-07-28T105355.438230.snakemake.log
```

``` sh
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% grep Error .snakemake/log/2023-07-28T105355.438230.snakemake.log  | awk '/executing/ {err[$4] = sprintf("%s%s", err[$4], $10)} END {for(i in err) print i, err[i]}'
genomics_db2vcf_scattered 2525063,2525064,
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--%
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--%
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--%
(snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% sacct -j 2525063,2525064
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode
------------ ---------- ---------- ---------- ---------- ---------- --------
2525063      smk-genom+     amilan csu-gener+          2  COMPLETED      0:0
2525063.bat+      batch            csu-gener+          2  COMPLETED      0:0
2525063.ext+     extern            csu-gener+          2  COMPLETED      0:0
2525064      smk-genom+     amilan csu-gener+          2  COMPLETED      0:0
2525064.bat+      batch            csu-gener+          2  COMPLETED      0:0
2525064.ext+     extern            csu-gener+          2  COMPLETED      0:0
```

I checked the log of one of those and GATK clearly reported that it had
finished.

When I did a dry-run it showed 500 or so jobs to go. But there are still
a few jobs running in SLURM:

    (snakemake-7.30.1) [login11: mega-non-model-wgs-snakeflow]--% myjobs 40
           JOBID PARTITION                                     NAME       USER ST            TIME  NODES   NODELIST(REASON)  CPUS   MIN_MEMORY   TIME_LIMIT    TIME_LEFT PRIORITY
         2524840    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:31:13      1    c3cpu-c11-u15-2     1        3600M   1-00:00:00     20:28:47 0.00000394601375
         2524820    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:31:51      1    c3cpu-c11-u13-1     1        3600M   1-00:00:00     20:28:09 0.00000394601375
         2524710    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:36:02      1      c3cpu-c9-u3-3     1        3600M   1-00:00:00     20:23:58 0.00000393949449
         2524716    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:36:02      1     c3cpu-c11-u1-1     1        3600M   1-00:00:00     20:23:58 0.00000393949449
         2524626    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:39:07      1      c3cpu-c9-u1-4     1        3600M   1-00:00:00     20:20:53 0.00000393926166
         2524629    amilan smk-make_gvcf_sections-bqsr_round=0,samp eriq@colos  R         3:39:07      1      c3cpu-c9-u1-4     1        3600M   1-00:00:00     20:20:53 0.00000393926166
         2525981    amilan smk-genomics_db_import_chromosomes-bqsr_ eriq@colos  R         1:58:18      1      c3cpu-c9-u1-1     2        7480M   1-00:00:00     22:01:42 0.00000393856317
         2525942    amilan smk-genomics_db_import_chromosomes-bqsr_ eriq@colos  R         2:04:23      1      c3cpu-c9-u1-2     2        7480M   1-00:00:00     21:55:37 0.00000393856317
         2525803    amilan smk-genomics_db_import_chromosomes-bqsr_ eriq@colos  R         2:19:09      1      c3cpu-c9-u1-1     2        7480M   1-00:00:00     21:40:51 0.00000393856317
         2527010    amilan smk-genomics_db_import_scaffold_groups-b eriq@colos  R         1:01:49      1      c3cpu-c9-u1-1     2       11000M   1-00:00:00     22:58:11 0.00000393716619
         2527016    amilan smk-genomics_db_import_chromosomes-bqsr_ eriq@colos  R         1:00:18      1      c3cpu-c9-u1-2     2        7480M   1-00:00:00     22:59:42 0.00000393693335
         2526759    amilan smk-genomics_db_import_chromosomes-bqsr_ eriq@colos  R         1:28:07      1      c3cpu-c9-u1-1     2        7480M   1-00:00:00     22:31:53 0.00000393693335

So, that is a weird error. I am going to let the currently running jobs
finish, and then check with sacct that they finished OK.
