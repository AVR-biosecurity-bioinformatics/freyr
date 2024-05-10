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

rt_fcid <- # convert Groovy to R list format
    stringr::str_extract_all(rt_fcid, pattern = "[^\\s,\\[\\]]+") %>% unlist()

### run R code

## join sample-level read tracking files into a single tibble
sample_tibble <- tibble() # new tibble
for (i in 1:length(rt_samples)) { # loop through .csv and add values to tibble as new rows
    new_csv <- read_csv(rt_samples[i], show_col_types = F, col_names = F)
    sample_tibble <- rbind(sample_tibble, new_csv)
}

colnames(sample_tibble) <- c("stage","sample_id","fcid","pcr_primers","fwd","rev") 

sample_tibble <- sample_tibble %>%
    mutate(sample_id = paste0(sample_id,"_",pcr_primers)) %>% # make sample_id consistent with "sample_id" + "pcr_primers" format
    dplyr::arrange(sample_id, pcr_primers, desc(fwd))

write_csv(sample_tibble, "sample_tibble.csv") # for debugging

## join fcid-level read tracking files into a single tibble

fcid_tibble <- tibble() # new tibble
for (i in 1:length(rt_fcid)) { # loop through .csv and add values to tibble as new rows
    new_csv <- read_csv(rt_fcid[i], show_col_types = F)
    fcid_tibble <- rbind(fcid_tibble, new_csv)
}

fcid_tibble <- fcid_tibble %>%
    mutate(sample_id = paste0(sample_id,"_",pcr_primers)) %>% # make sample_id consistent with "sample_id" + "pcr_primers" format
    dplyr::arrange(sample_id, pcr_primers, desc(pairs))

write_csv(fcid_tibble, "fcid_tibble.csv") # for debugging
