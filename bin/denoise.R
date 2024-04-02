#!/usr/bin/env Rscript

# check variables defined


### run R code
if (direction == "forward") { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if (direction == "reverse") {
    direction_short <- "R"
} else {
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

if (n_pass == "first" && priors == "null") { # first pass condition; no priors
    priors <- NA

    # run dada2
    dada2::dada(
    derep = reads, 
    err = errormodel, 
    multithread = FALSE, 
    priors = priors, 
    selfConsist = FALSE, 
    pool = FALSE, 
    verbose = TRUE
)
} else if (n_pass == "second" && !priors == "null") { # second pass condition; priors included

    # run dada2
    dada2::dada(
    derep = reads, 
    err = errormodel, 
    multithread = FALSE, 
    priors = priors, 
    selfConsist = FALSE, 
    pool = FALSE, 
    verbose = TRUE
)
} else {
    stop(" 'n_pass' variable must be 'first' or 'second', and priors must be 'null' or defined! ")
}

# alex function
# step_dada2_single2(
#     fcid, 
#     sample_id, 
#     input_dir,
#     pcr_primers, 
#     output, 
#     qc_dir, 
#     error_model, 
#     read="F",
#     priors = NA, 
#     quiet=FALSE, 
#     multithread=FALSE
# )

