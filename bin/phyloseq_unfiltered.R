#!/usr/bin/env Rscript


## check and define variables
merged_tax_list <- # convert Groovy to R list format
    stringr::str_extract_all(merged_tax, pattern = "[^\\s,\\[\\]]+") %>% unlist()

merged_tax_list <- lapply(merged_tax_list, readRDS) # read in taxtabs and store as list of matrices

### run R code
taxtab <- merged_tax_list %>% 
    purrr::map(~{ .x %>% tibble::as_tibble(rownames = "OTU") }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct() # Remove any exact duplicates (unlikely as different primers used, but possible)

write_csv(taxtab, "taxtab.csv")


# to run classic pipeline code, need:
# - merged seqtab (all sequences across all flowcells and loci)
# - merged taxtable (all sequences across all flowcells and loci, identified)
# - samplesheet (ie. samdf)
# - loci params csv file (to be used inside this process)

# should input all the files independently, merge within this code, then run filtering etc.
# tax tables (output of MERGE_TAX) are already per locus (merged across flowcells)
# sequence tables (output of FILTER_SEQTAB) are per locus x flowcell

## merge sequence tables across flowcells and loci



## merge taxonomy tables across loci

## 

## runs step_phyloseq() on merged seqtabs, taxtabs and samplesheet

## runs step_output_summary() and step_output_ps() on unfiltered data

## runs rareplot() and saves plot

## runs step_filter_phyloseq() on each locus independently

## runs merge_phyloseq_new() from functions.R

## runs step_output_summary() and step_output_ps() on filtered data





# stop(" *** stopped manually *** ") ##########################################
