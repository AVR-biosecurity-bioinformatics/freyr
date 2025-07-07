#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

tax_summary_list            <- args$tax_summary_list

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "readr",
    "stringr",
    "tibble",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

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

readr::write_csv(ts_merged,"taxonomic_assignment_summary_combined.csv") # write to .csv

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})