#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dada2",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "direction",
    "fcid",
    "pcr_primers",
    "sample_id",
    "reads",
    "errormodel",
    "n_pass",
    "priors"
)
lapply(nf_vars, nf_var_check)

### run R code
errormodel <- readRDS(errormodel) # import error model

if ( direction == "forward" ) { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if ( direction == "reverse" ) {
    direction_short <- "R"
} else if ( direction == "single" ) {
    direction_short <- "S"
} else {
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

if ( n_pass == "first" && priors == "NO_FILE" ) { # first pass condition; no priors
    priors <- character(0)

    # run dada2
    dada_output <- dada2::dada(
        derep = reads, 
        err = errormodel, 
        multithread = FALSE, 
        priors = priors, 
        selfConsist = FALSE, 
        pool = FALSE, 
        verbose = TRUE
    )

    saveRDS(dada_output, paste0(sample_id,"_",pcr_primers,"_dada1",direction_short,".rds"))

} else if ( n_pass == "second" && !priors == "NO_FILE" ) { # second pass condition; priors included
    # import priors
    priors <- readRDS(priors)

    # run dada2
    dada_output <- dada2::dada(
        derep = reads, 
        err = errormodel, 
        multithread = FALSE, 
        priors = priors, 
        selfConsist = FALSE, 
        pool = FALSE, 
        verbose = TRUE
    )

    saveRDS(dada_output, paste0(sample_id,"_",pcr_primers,"_dada2",direction_short,".rds"))

} else {
    stop(" 'n_pass' variable must be 'first' or 'second', and priors must be 'NO_FILE' or defined! ")
}

