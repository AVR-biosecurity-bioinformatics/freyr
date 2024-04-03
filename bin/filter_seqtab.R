#!/usr/bin/env Rscript

# check variables defined

### run R code
mergers_list <- # convert input mergers list from Groovy format to R format
    stringr::str_extract_all(
        mergers, 
        pattern = "\\S+?\\.rds" 
        ) %>% 
    unlist()

print(mergers_list)

# stop(" *** stopped manually *** ") ##########################################

# Construct sequence table for fcid x pcr_primers from merged reads per sample
seqtab <- dada2::makeSequenceTable(mergers)
saveRDS(seqtab, paste0(sample_id,"_",pcr_primers,"_seqtab.rds"))

