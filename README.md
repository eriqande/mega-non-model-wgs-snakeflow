mega-non-model-wgs-snakeflow
================

## Quick install and run

If you would like to put this on your system and test it running on a
single node (more later about using SLURM for deployment across multiple
nodes) you have to clone this repository and then download the
pseudo-genome used for the included test data set (in `.test`).

You must have Snakemake (version &gt; 6.0) in the active environment.

In short, here are the steps to install and run the `.test`.

``` sh
# clone the repo
git clone git@github.com:eriqande/mega-non-model-wgs-snakeflow.git

# download the tarball with the genome in it and then move that
# into resources/
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1LMK-DCkH1RKFAWTR2OKEJ_K9VOjJIZ1b' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1LMK-DCkH1RKFAWTR2OKEJ_K9VOjJIZ1b" -O non-model-wgs-example-data.tar && rm -rf /tmp/cookies.txt

# untar the tarball
tar -xvf non-model-wgs-example-data.tar

# move the genome from the extracted tarball into mega-non-model-wgs-snakeflow/resources/

cp non-model-wgs-example-data/resources/genome.fasta mega-non-model-wgs-snakeflow/resources/
```

## What the user must do and values to be set, etc

-   Choose an Illuminaclip adapter fasta (in config)

## Assumptions

-   Paired end

## Things fixed or added relative to JK’s snakemake workflow

-   fastqc on both reads
-   don’t bother with single end
-   add adapters so illumina clip can work
-   benchmark each rule
-   use genomicsDBimport
-   allow for merging of lots of small scaffolds into genomicsDB
