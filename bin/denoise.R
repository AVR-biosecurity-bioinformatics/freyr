#!/usr/bin/env Rscript

# check variables defined

### run R code
errormodel <- readRDS(errormodel) # import error model

if ( direction == "forward" ) { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if ( direction == "reverse" ) {
    direction_short <- "R"
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

