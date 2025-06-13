#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "phyloseq",
    "purrr",
    "readr",
    "scales",
    "stringr",
    "tibble",
    "data.table",
    "tidyr",
    "vegan",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "primers",
    "ps_file",
    "min_sample_reads"
)
lapply(nf_vars, nf_var_check)

## check and define variables
ps <- readRDS(ps_file)
min_sample_reads <- as.numeric(min_sample_reads)

### run R code

## creates accumulation curve plots and saves plot
gg.acc_curve <- 
    rareplot(
        ps, 
        step=1L, 
        threshold = min_sample_reads
    )

pdf(file=paste0("accumulation_curve_",primers,".pdf"), width = 11, height = 8 , paper="a4r")
    print(gg.acc_curve)
try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################
