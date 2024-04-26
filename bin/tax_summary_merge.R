#!/usr/bin/env Rscript


## check and define variables
print(tax_summary_list)
tax_summary_list <- # convert Groovy to R list format
    stringr::str_extract_all(tax_summary_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()
print(tax_summary_list)
tax_summary_list <- lapply(tax_summary_list, readRDS) # read in taxtabs and store as list of tibbles

print(tax_summary_list)

### run R code


stop(" *** stopped manually *** ") ##########################################

