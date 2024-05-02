#!/usr/bin/env Rscript


## check and define variables
taxtab <- readRDS(taxtab)

seqtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(seqtab_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()
seqtab_list <- lapply(seqtab_list, readRDS) # read in seqtabs and store as list of matrices

samdf <- readr::read_csv(samdf, show_col_types = FALSE)

### run R code

# to run classic pipeline code, need:
# - merged seqtab (all sequences across all flowcells and loci)
# - merged taxtable (all sequences across all flowcells and loci, identified)
# - samplesheet (ie. samdf)
# - loci params csv file (to be used inside this process)

# should input all the files independently, merge within this code, then run filtering etc.
# tax tables (output of MERGE_TAX) are already per locus (merged across flowcells)
# sequence tables (output of FILTER_SEQTAB) are per locus x flowcell

# ## merge list of taxtabs
# taxtab <- taxtab_list %>% 
#     purrr::map(~{ .x %>% tibble::as_tibble(rownames = "OTU") }) %>% # convert each matrix to tibble with rownames to OTU column
#     dplyr::bind_rows() %>% # bind tibbles into one
#     dplyr::distinct() # remove any exact duplicate rows (unlikely as different primers used per input tibble)

taxtab %>% # save for debugging
    tibble::as_tibble(rownames = "OTU") %>%
    write_csv(., paste0("taxtab_", pcr_primers, ".csv")) 

## merge sequence tables across flowcells and loci
if ( length(seqtab_list) > 1 ){ # if there is more than one seqtab, merge together
    seqtab_final <- dada2::mergeSequenceTables(tables=seqtab_list)
} else if( length(seqtab_list) == 1 ) { # if there is only one seqtab, keep and unlist
    seqtab_final <- seqtab_list %>% unlist()
}

seqtab_final %>% # save for debugging
    tibble::as_tibble(rownames = "OTU") %>% 
    write_csv(., paste0("seqtab_final_", pcr_primers, ".csv"))

## runs step_phyloseq() on merged seqtabs, taxtabs and samplesheet

ps <-   step_phyloseq(
            seqtab = seqtab_final,
            taxtab = taxtab,
            samdf = samdf,
            seqs = NULL,
            phylo = NULL,
            name_variants = FALSE
        )

## name OTUs using hash
taxa_names(ps) <- tax_table(ps)[,ncol(tax_table(ps))]
tax_table(ps) <- tax_table(ps)[,1:ncol(tax_table(ps))-1] # remove hash 'rank' from taxtab


saveRDS(ps, paste0("ps_",pcr_primers,".rds"))




# stop(" *** stopped manually *** ") ##########################################

## runs step_output_summary() and step_output_ps() on unfiltered data

## runs rareplot() and saves plot

## runs step_filter_phyloseq() on each locus independently

## runs merge_phyloseq_new() from functions.R

## runs step_output_summary() and step_output_ps() on filtered data





# stop(" *** stopped manually *** ") ##########################################
