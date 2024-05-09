#!/usr/bin/env Rscript

# code produces:
# - plot of N reads retained after each step, faceted by flowcell, stacked bars coloured by locus
# - .csv file that show same info by per sample (split by locus)

### required inputs:
# - unfiltered phyloseq
# - filtered phyloseq
# - reads input for trimming
# - 


## check and define variables 
rt_samples <- # convert Groovy to R list format
    stringr::str_extract_all(rt_samples, pattern = "[^\\s,\\[\\]]+") %>% unlist()

samples_tibble <- tibble() # new tibble
for (i in 1:length(rt_samples)) { # loop through .rds files, adding distinct sequences to tibble
    new_csv <- read_csv(rt_samples, show_col_types = F, col_names = F)
    samples_tibble <- rbind(samples_tibble, new_csv)
}

write_csv(samples_tibble, "samples_tibble.csv")

### run R code
