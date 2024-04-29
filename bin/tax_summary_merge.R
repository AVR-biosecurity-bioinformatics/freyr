#!/usr/bin/env Rscript


## check and define variables

tax_summary_list <- # convert Groovy to R list format
    stringr::str_extract_all(tax_summary_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()

tax_summary_list <- lapply(tax_summary_list, readRDS) # read in taxtabs and store as list of tibbles



### run R code
summary_full <- tibble() # create empty tibble
for (i in 1:length(tax_summary_list)){ # loop through list of tibbles
    summary_full <- bind_rows(summary_full, tax_summary_list[i]) # combine tibbles into one tibble
}

summary_full <- summary_full %>%
    distinct() %>%  # remove identical rows, leaving one (ie. duplicate sequences between flow cells)
    arrange(pcr_primers, OTU_hash) # sort final tibble by pcr_primer and OTU_hash

write_csv(summary_full,"taxonomic_assignment_summary.csv") # write to .csv

# stop(" *** stopped manually *** ") ##########################################

