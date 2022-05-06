# This rips through the fai file and creates two files that
# we can then use to break up the genotyping over chromosomes
# and also over merged blobs of shorter scaffolds.



#save(snakemake, file = "/tmp/smk.rda")
#load(file = "/tmp/smk.rda")

# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# this should put the search paths in the correct order, otherwise
# the locally installed libraries take precendence with disastrous
# consequences for reproducibility
LP <- .libPaths()
LP <- LP[grep("conda", LP)]
.libPaths(new = LP)

library(tidyverse)


fai_file <- snakemake@input$fai
chromo_file <- snakemake@output$chrom
scaffo_file <- snakemake@output$scaff

# chromosome-prefix
CP <- snakemake@config$ref$chromosome_prefix
scaff_frac <- as.numeric(snakemake@config$ref$scaff_frac)

fai <- read_tsv(
  fai_file, 
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
  write_tsv(chromo_file)

# now, figure out the mean length of the chromosomes.
# We will set our scaff groups to cumulatively be scaff_frac of that.
mean_cl <- chromos %>% pull(len) %>% mean()

bin_length <- as.integer(mean_cl * scaff_frac)

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
  write_tsv(scaffo_file)



