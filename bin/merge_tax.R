#!/usr/bin/env Rscript


## check and define variables

taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(taxtab, pattern = "[^\\s,\\[\\]]+") %>% unlist()

taxtab_list <- lapply(taxtab_list, readRDS) # read in taxtabs and store as list of tibbles

### run R code

tax_merged <- taxtab_list %>% 
    purrr::map(~{ .x %>% tibble::as_tibble(rownames = "OTU") }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() # Remove any exact duplicates being in different seqtab

# create hash of OTU sequence that can be used as OTU name
# 'OTU_hash' column is now the 'name' of the sequence
OTU_hash <- lapply(tax_merged$OTU, rlang::hash) %>% unlist()
tax_merged <- cbind(tax_merged, OTU_hash)

# saveRDS(tax_merged, paste0(pcr_primers,"_tax_merged.rds"))

# Check for duplicated ASVs across taxtabs
if(any(duplicated(tax_merged$OTU))){
    warning("Duplicated ASVs detected, selecting first occurrence")
    merged_tax <- tax_merged %>%
        dplyr::group_by(OTU) %>%
        dplyr::slice(1) %>%
        tibble::column_to_rownames("OTU") %>%
        as.matrix()
} else{
    merged_tax <- tax_merged %>%
    tibble::column_to_rownames("OTU") %>%
    as.matrix()
}

# print(out)

### TODO: Not sure dimensions of merged seqtab
# Check that output dimensions match input
### .y is mergedseqtab
# if(!all(rownames(out) %in% colnames(.y))){
#     stop("Number of ASVs classified does not match the number of input ASVs")
# }

saveRDS(merged_tax, paste0(pcr_primers,"_merged_tax.rds"))

# stop(" *** stopped manually *** ") ##########################################
