#!/usr/bin/env Rscript


## check and define variables

tax_summary_list <- # convert Groovy to R list format
    stringr::str_extract_all(tax_summary_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()

tax_summary_list <- lapply(tax_summary_list, readRDS) # read in taxtabs and store as list of tibbles



### run R code
summary_full <- tibble()
for (i in 1:length(tax_summary_list)){
    summary_full <- bind_rows(summary_full, tax_summary_list[i])
}

summary_full <- summary_full %>%
    distinct() 

write_csv(summary_full,"taxonomic_assignment_summary.csv")

# stop(" *** stopped manually *** ") ##########################################

