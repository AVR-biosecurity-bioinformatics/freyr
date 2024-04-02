#!/usr/bin/env Rscript

# check variables defined

### run R code
priors_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        priors, 
        pattern = "\\S+?\\.rds" 
        ) %>% 
    unlist()

priors_seq <- list()

for (i in 1:length(priors_list)) {
    seq_tmp <- readRDS(priors_list[i])$sequence
    priors_seq <- append(priors_seq, list(seq_tmp))
}

as.data.frame(priors_seq)

print(priors_seq)


## keep only duplicates within each component

## keep only duplicates between components (ie. unnest then check dupes)

stop(" *** stopped manually *** ") ##########################################

# Only keep the ones that appear across more than one sample
priors <- unlist(process$priors)
priors <- names(table(priors))[table(priors) > 1]