

library(tidyverse)


times <- list(
  scratch = read_tsv("development/scratch-timings-etc/scratch-benchmarks.bmk"),
  nfs = read_tsv("development/scratch-timings-etc/nfs-benchmarks.bmk")
) %>%
  bind_rows(.id = "file_system") %>%
  select(-`h:m:s`) %>%
  rename(wall_time_seconds = s) %>%
  pivot_longer(cols = wall_time_seconds:cpu_time) %>%
  separate(file, into = c("rule", "sample"), sep = "/") %>%
  pivot_wider(
    names_from = file_system,
    values_from = value
  )

plot_it <- function(X) {
  X %>%
  ggplot(aes(x = nfs, y = scratch, fill = rule)) +
  geom_point(shape = 21, stroke = 0.1) + 
  facet_grid(rule ~ name, scales = "free_x") +
  geom_abline(slope = 1, intercept = 0)
}

# now plot those, facet-gridding on rule and name
times %>%
  filter(name == "wall_time_seconds", rule %in% c("fastqc_read1", "fastqc_read2")) %>%
  plot_it()
