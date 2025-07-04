#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
taxtabs                     <- args$taxtabs

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "purrr",
    "readr",
    "rlang",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables
taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(taxtabs, pattern = "[^\\s,\\[\\]]+") %>% unlist()

taxtab_list <- lapply(taxtab_list, readr::read_csv) # read in taxtabs and store as list of tibbles

### run R code

tax_merged <- 
    taxtab_list %>% 
    dplyr::bind_rows() %>%
    dplyr::distinct() # Remove any exact duplicates being in different tables

# Check for duplicated ASVs across taxtabs
if(any(duplicated(tax_merged$seq_name))){
    warning("Duplicated ASVs detected, selecting first occurrence")
    merged_tax <- 
        tax_merged %>%
        dplyr::group_by(seq_name) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
} else{
    merged_tax <- tax_merged
}


readr::write_csv(merged_tax, paste0(primers,"_merged_tax.csv"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})