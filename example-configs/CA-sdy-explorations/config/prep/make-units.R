

library(tidyverse)

paths <- read_tsv("example-configs/CA-sdy-explorations/config/prep/paths.tsv")
meta <- read_csv("example-configs/CA-sdy-explorations/config/prep/wgs-chinook-samples.csv") %>%
  select(vcf_name, NMFS_DNA_ID) %>%
  rename(sample_id = NMFS_DNA_ID)


units <- paths %>%
  extract(
    path, 
    into = c("vcf_name", "read"),
    regex = ".*(DPCh_plate[0-9]_.+_.+)_R([12]).*gz",
    remove = FALSE,
    convert = TRUE
  ) %>%
  rename(fq = path) %>%
  pivot_wider(names_from = read, values_from = c(fq, kb), names_sep = "") %>%
  mutate(
    sample = sprintf("s%03d", 1:n()),
    unit = 1,
    library = "lib",
    flowcell = "whatever",
    platform = "ILLUMINA",
    lane = 1,
    barcode = vcf_name
  ) %>%
  left_join(meta , by = "vcf_name") %>%
  select(sample, unit, library, flowcell, platform, barcode, sample_id, lane, fq1, fq2, kb1, kb2)


write_tsv(units, file = "example-configs/CA-sdy-explorations/config/units.tsv")
