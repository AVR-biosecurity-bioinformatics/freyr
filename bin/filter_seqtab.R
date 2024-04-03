#!/usr/bin/env Rscript

# check variables defined

### run R code


# Construct sequence table for fcid x pcr_primers from merged reads per sample
seqtab <- dada2::makeSequenceTable(mergers)
saveRDS(seqtab, paste0(sample_id,"_",pcr_primers,"_seqtab.rds"))

