#!/usr/bin/env Rscript

# check variables defined


### run R code
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

if (high_sen_mode == "yes" && !priors == "null") {

} else if (high_sen_mode == "no" && priors == "null") {

} else {
    stop(" 'high_sen_mode' variable must be 'yes' or 'no'! ")
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

# direct function from dada2
dada2::dada(
    derep = reads, 
    err = errormodel, 
    multithread = FALSE, 
    priors = priors, 
    selfConsist = FALSE, 
    pool = FALSE, 
    verbose = TRUE
)