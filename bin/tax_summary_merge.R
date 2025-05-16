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
# read in taxtabs and store as list of tibbles
ts_merged <- 
    tax_summary_list %>% 
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(seq_name) %>%
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(seq_name)

readr::write_csv(ts_merged,"taxonomic_assignment_summary.csv") # write to .csv

# stop(" *** stopped manually *** ") ##########################################
