# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)


# these are while developing
#chroms_file = ".test/config/chromosomes.tsv"
#scaffs_file = ".test/config/scaffold_groups.tsv"
#max_chunk <- 900000  # shooting for just < than 1 Mb chunks


chroms_file <- snakemake@input$chroms
scaffs_file <- snakemake@input$scaffs
max_chunk <- as.integer(snakemake@params$binsize)
outfile <- snakemake@output$interval_list


# read the files in and prepare them so they have
# columns that are comparable
chroms_raw <- read_tsv(chroms_file) 
scaffs_raw <- read_tsv(scaffs_file)

# only do the following if the chromosomes.tsv file is more than just the header
if(nrow(chroms_raw) > 0) {
  chroms <- chroms_raw %>%
    mutate(id = chrom, .before = chrom) %>%
    group_by(id) %>%
    mutate(cumul = cumsum(num_bases)) %>%
    ungroup()
  
  
  
  
  # Now, we handle the chromosomes and the scaffold groups a little
  # differently.  (It's a shame really, they should just all be scaff_groups,
  # with the chromosomes just being scaff groups of size 1.  Sigh...)
  
  # for the chromos, we just find a good length to chop them up into.  
  # Here is a function that tries to get them all close to an equal length
  # that is less than max_chunk:
  #' @param L length of the chromosome
  #' @param m max length of chunk
  #' @param S the starting value
  #' @return returns a tibble with start and stop columns of the different
  #' chunks. (base-1 indexed, inclusive).
  chromo_chopper = function(L, m, S = 1) {
    # The number of bins will be:
    B = ceiling(L / m)
    
    # if m divides L perfectly, mstar is m and last_one = mstar
    if((L %% m) == 0) {
      mstar = m
      last_one = mstar
    } else {  # otherwise
      mstar = ceiling(L / B)
      last_one = L - mstar * (B - 1)
    }
    
    starts <- seq(S, S + B * (mstar - 1), by = mstar)
    ends <- starts + mstar - 1
    ends[length(ends)] <- starts[length(ends)] + last_one - 1
    
    tibble(
      start = starts,
      end = ends
    ) %>%
      mutate(scatter_length = end - start + 1) %>%
      mutate(scatter_idx = sprintf("scat_%04d", 1:n()), .before = start)
  }
  
  
  chroms2 <- chroms %>%
    group_by(id) %>%
    mutate(scats = map(num_bases, chromo_chopper, m = max_chunk)) %>%
    unnest(scats) %>%
    select(id, scatter_idx, chrom, start, end, scatter_length)
} 



# Only to the following if the scaff_groups file is more than just the header.
if(nrow(scaffs_raw) > 0) {
  
  scaffs <- scaffs_raw  %>%
    rename(num_bases = len)
  
  # now, for the scaffold groups we will do it a little differently.
  # Any scaffold that is greater than the max_chunk will get broken down into
  # roughly equal-sized chunks that are all less than the max_chunk.
  scaff_chopper <- function(L, m) {
    if(L <= m) return(
      tibble(
        start = 1,
        end = L,
        seg_length = L
      )
    )
    else {
      return(
        chromo_chopper(L, m) %>%
          rename(seg_length = scatter_length) %>%
          select(-scatter_idx)
      )
    }
  }
  
  scaffs2 <- scaffs %>%
    ungroup() %>%
    mutate(chops = map(num_bases, scaff_chopper, m = max_chunk)) %>%
    unnest(chops) %>%
    group_by(id) %>%
    mutate(cumul = cumsum(seg_length))
  
  
  # now we will group by id and within each one, just iterate through
  # the number of bases and assign to different scat groups.
  #' @param s the number of bases in each segment (i.e. the seg_length)
  #' @param m the desired max chunk size
  iteratively_assign_scats <- function(n, m) {
    c <- 0  # cumulative in scat group
    scg <- 1
    ret <- character(length(n))
    for(i in 1:length(n)) {
      c <- c + n[i]
      if(c > m) {
        scg <- scg + 1
        c <- n[i]
      }
      ret[i] <- sprintf("scat_%04d", scg)
    }
    ret
  }
  
  scaffs3 <- scaffs2 %>%
    group_by(id) %>%
    mutate(scatter_idx = iteratively_assign_scats(n = seg_length, m = max_chunk)) %>%
    group_by(id, scatter_idx) %>%
    mutate(scatter_length = sum(seg_length)) %>%
    ungroup() %>%
    select(id, scatter_idx, chrom, start, end, scatter_length)
}

if(nrow(scaffs_raw) > 0 && nrow(chroms_raw) > 0) {
  chrom_and_sg <- bind_rows(chroms2, scaffs3)
} else if(nrow(chroms_raw) > 0) {
  chrom_and_sg <- chroms2
} else if(nrow(scaffs_raw) > 0) {
  chrom_and_sg <- scaffs3
} else {
  stop("Error. At least one of chromosomes or scaffold_groups must have entries.")
}

write_tsv(chrom_and_sg, file = snakemake@output$tsv)
