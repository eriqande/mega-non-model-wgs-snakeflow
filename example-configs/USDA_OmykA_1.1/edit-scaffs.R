
library(tidyverse)

# just some stuff I am doing to get things good to go
sg <- read_tsv("example-configs/non-RNAsed-trout/scaffold_groups.tsv")


# get things so that we don't have more than 30 in each scaffold group
sg_new <- sg %>%
  mutate(
    num = 1:n(),
    gp = floor((num - 1) / 30) + 1 
  ) %>%
  mutate(id = sprintf("scaff_group%03d", gp)) 


write_tsv(sg_new, file = "example-configs/USDA_OmykA_1.1/scaffold_groups.tsv")
