#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    # "Biostrings",
    # "bs4Dash",
    # "clustermq",
    "dada2",
    # "DECIPHER",
    # "dplyr",
    # "future",
    "ggplot2",
    # "gridExtra",
    # "gt",
    # "magrittr",
    # "markdown",
    # "ngsReports",
    # "patchwork",
    # "phyloseq",
    # "pingr",
    # "purrr",
    "readr",
    # "rlang",
    # "rstudioapi",
    # "savR",
    # "scales",
    # "seqateurs",
    # "shiny",
    # "shinybusy",
    # "shinyWidgets",
    # "ShortRead",
    "stringr",
    # "taxreturn",
    # "tibble",
    # "tidyr",
    # "vegan",
    # "visNetwork",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

# check variables defined


### run R code (from step_errormodel)
if (!direction %in% c("forward","reverse")) { 
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

if (direction == "forward") {
    direction_short <- "F"
} else if (direction == "reverse") {
    direction_short <- "R"
} else {
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

reads_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()

# input_dir <- normalizePath(input_dir)
# output <- normalizePath(output)
# qc_dir <- normalizePath(qc_dir)

if(direction == "forward"){
    # filts <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
    message(paste0("Modelling forward read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else if (direction == "reverse"){
    # filts <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
    message(paste0("Modelling reverse read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
} else {
    stop ("Read direction must be 'forward' or 'reverse'!")
}

## Learn error rates 
err <- dada2::learnErrors(
        reads_list, 
        multithread = FALSE, 
        nbases = 1e8,
        randomize = FALSE, 
        qualityType = "FastqQuality",
        verbose=TRUE
        )
    
### TODO: increase nbases or use default (1e8)

## save output as .rds file
saveRDS(err, paste0(fcid,"_",pcr_primers,"_errormodel",direction_short,".rds"))

## write out errors for diagnostics
readr::write_csv(as.data.frame(err$trans), paste0(fcid,"_",pcr_primers,"_err",direction_short,"_observed_transitions.csv"))
readr::write_csv(as.data.frame(err$err_out), paste0(fcid,"_",pcr_primers,"_err",direction_short,"_inferred_errors.csv"))

## output error plots to see how well the algorithm modelled the errors in the different runs
p1 <- dada2::plotErrors(err, nominalQ = TRUE) + ggtitle(paste0(pcr_primers, " ", fcid," ",direction," Reads"))
pdf(paste0(fcid,"_", pcr_primers, "_", direction_short,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
    plot(p1)
try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################