#!/usr/bin/env Rscript


## check and define variables

#### INCOMPLETE
fcid_list <- # convert Groovy to R list format
    stringr::str_extract_all(fcid, pattern = "[^\\s,\\[\\]]+") %>% unlist()

meta_list <-  # convert Groovy to R list format
    stringr::str_extract_all(meta, pattern = "[^\\[\\]]+") %>% 
    unlist() %>% str_subset(pattern = "^[^,]") # remove elements that start with ","
meta_list2 <- meta_list %>% stringr::str_split(pattern = ", ") # nested list; elements are samples x pcr_primers
#### INCOMPLETE above


taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(taxtab, pattern = "[^\\s,\\[\\]]+") %>% unlist()

taxtab_list <- lapply(taxtab_list, readRDS) # read in taxtabs and store as list of tibbles

saveRDS(taxtab_list, paste0(pcr_primers,"_taxtab_list.rds"))

### run R code

tax_merged <- taxtab_list %>% 
    purrr::map(~{ .x %>% tibble::as_tibble(rownames = "OTU") }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() # Remove any exact duplicates from save ASV being in different seqtab

saveRDS(tax_merged, paste0(pcr_primers,"_tax_merged.rds"))

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
