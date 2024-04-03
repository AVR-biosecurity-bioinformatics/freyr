#!/usr/bin/env Rscript

# check variables defined

### run R code
seqtab_list <- # convert input seqtab list from Groovy format to R format
    stringr::str_extract_all(
        seqtab, 
        pattern = "\\S+?\\.rds" 
        ) %>% 
    unlist()

print(seqtab_list)

merged_seqtab <- 
    dada2::mergeSequenceTables(
        tables = seqtab_list,
        repeats = "error",
        orderBy = "abundance",
        tryRC = FALSE
        )

print(merged_seqtab)

stop(" *** stopped manually *** ") ##########################################


