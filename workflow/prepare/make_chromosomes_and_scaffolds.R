

# This rips through the fai file and creates two files that
# we can then use to break up the genotyping over chromosomes
# and also over merged blobs of shorter scaffolds.

# This is intended to be done interactively by the user, and will
# create the files chromosomes.tsv and scaffold_groups.tsv

library(tidyverse)

# Prefix that the true chromosomes start with.  All others are assumed
# to be "scaffolds" and will be merged.

# Note that I use this prefix matching approach to define "chromosomes" and
# "scaffolds", but you could use any other means of partitioning them. All that
# matters is that you have mutually exclusive sets in the variables
# chromos and scaffs down below
CP <- "CM0"

fai <- read_table(
  "resources/genome.fasta.fai", 
  col_names = c("chrom", "len", "x", "y", "z")
) %>%
  select(chrom, len)
  

chromos <- fai %>%
  filter(str_detect(chrom, str_c("^", CP)))

scaffs <- fai %>%
  filter(!str_detect(chrom, str_c("^", CP)))

# now, write chromosomes.tsv
chromos %>%
  mutate(id = sprintf("chromo%03d", 1:n()), .before = chrom) %>%
  select(-len) %>%
  write_tsv("chromosomes.tsv")

# now, figure out the mean length of the chromosomes.
# We will set our scaff groups to cumulatively be 40% of that.
mean_cl <- chromos %>% pull(len) %>% mean()

bin_length <- as.integer(mean_cl * 0.40)

scaffs %>%
  mutate(
    cumul = cumsum(len),
    part = as.integer(cumul / bin_length) + 1L
  ) %>%
  mutate(
    id = sprintf("scaff_group%03d", part),
    .before = chrom
  ) %>%
  select(-part) %>%
  write_tsv("scaffold_groups.tsv")
