#!/usr/bin/env Rscript

# check variables defined


### run R code

# alex function
step_dada2_single2(
    fcid, 
    sample_id, 
    input_dir,
    pcr_primers, 
    output, 
    qc_dir, 
    error_model, 
    read="F",
    priors = NA, 
    quiet=FALSE, 
    multithread=FALSE
)

# direct function from dada2
dada2::dada(
    filts, 
    err = error_model, 
    multithread = FALSE, 
    priors = priors, 
    selfConsist = FALSE, 
    pool = FALSE, 
    verbose = TRUE
)