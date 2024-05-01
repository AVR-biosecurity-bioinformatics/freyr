#!/usr/bin/env Rscript


## check and define variables
taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(taxtab_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()
taxtab_list <- lapply(taxtab_list, readRDS) # read in taxtabs and store as list of matrices

seqtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(seqtab_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()
seqtab_list <- lapply(seqtab_list, readRDS) # read in seqtabs and store as list of matrices


### run R code

# to run classic pipeline code, need:
# - merged seqtab (all sequences across all flowcells and loci)
# - merged taxtable (all sequences across all flowcells and loci, identified)
# - samplesheet (ie. samdf)
# - loci params csv file (to be used inside this process)

# should input all the files independently, merge within this code, then run filtering etc.
# tax tables (output of MERGE_TAX) are already per locus (merged across flowcells)
# sequence tables (output of FILTER_SEQTAB) are per locus x flowcell

## merge taxtables across loci
taxtab <- taxtab_list %>% 
    purrr::map(~{ .x %>% tibble::as_tibble(rownames = "OTU") }) %>% # convert each matrix to tibble with rownames to OTU column
    dplyr::bind_rows() %>% # bind tibbles into one
    dplyr::distinct() # remove any exact duplicate rows (unlikely as different primers used per input tibble)

write_csv(taxtab, "taxtab.csv") # write intermediate output table

taxtab_matrix <- taxtab %>% # convert tibble to matrix format
    tibble::column_to_rownames("OTU") %>%
    as.matrix()

## merge sequence tables across flowcells and loci
if ( length(seqtab_list) > 1){ # if there is more than one seqtab, merge together
    seqtab_final <- dada2::mergeSequenceTables(tables=seqtab_list)
} else if( length(seqtab_list) == 1) { # if there is only one seqtab, keep and unlist
    seqtab_final <- seqtab_list %>% unlist()
}

seqtab_final %>% # save for debugging
    tibble::as_tibble(rownames = "OTU") %>% 
    write_csv(., "seqtab_final.csv")




## runs step_phyloseq() on merged seqtabs, taxtabs and samplesheet

# step_phyloseq(
#     seqtab = seqtab_final,
#     taxtab = taxtab_matrix,
#     samdf = ,
#     seqs = NULL,
#     phylo = NULL,
#     name_variants = FALSE
# )

## runs step_output_summary() and step_output_ps() on unfiltered data

## runs rareplot() and saves plot

## runs step_filter_phyloseq() on each locus independently

## runs merge_phyloseq_new() from functions.R

## runs step_output_summary() and step_output_ps() on filtered data





# stop(" *** stopped manually *** ") ##########################################
