#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "readr",
    "stringr",
    "tibble",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "tax_summary_list"
)
lapply(nf_vars, nf_var_check)

## check and define variables

tax_summary_list <- # convert Groovy to R list format
    stringr::str_extract_all(tax_summary_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()

tax_summary_list <- lapply(tax_summary_list, readRDS) # read in taxtabs and store as list of tibbles


### run R code
summary_full <- tibble::tibble() # create empty tibble
for (i in 1:length(tax_summary_list)){ # loop through list of tibbles
    summary_full <- dplyr::bind_rows(summary_full, tax_summary_list[i]) # combine tibbles into one tibble
}

summary_full <- summary_full %>%
    dplyr::distinct() %>%  # remove identical rows, leaving one (ie. duplicate sequences between flow cells)
    dplyr::arrange(OTU_hash) # sort final tibble by OTU_hash

readr::write_csv(summary_full,"taxonomic_assignment_summary.csv") # write to .csv

# stop(" *** stopped manually *** ") ##########################################
