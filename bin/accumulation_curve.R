#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
ps_file                     <- args$ps_file
min_sample_reads            <- args$min_sample_reads

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

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

## check and define variables
ps <- readRDS(ps_file)
min_sample_reads <- as.numeric(min_sample_reads)

### run R code

## creates accumulation curve plots and saves plot
gg.acc_curve <- 
    rareplot(
        ps, 
        step="auto", 
        threshold = min_sample_reads
    )

pdf(file=paste0("accumulation_curve_",primers,".pdf"), width = 11, height = 8 , paper="a4r")
    print(gg.acc_curve)
try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})