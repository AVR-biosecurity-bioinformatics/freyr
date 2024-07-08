#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "ggplot2",
    "magrittr",
    "purrr",
    "readr",
    "savR",
    "seqateurs",
    "stringr",
    "tibble",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

#### TODO: Pull stats from top of index_switching.pdf and print to a final run report
### also flag index combos that have higher than average + (some stat)

#### TODO: use index_switch_calc.txt in jack_notes to run same process but in bash, which should be much faster

# check variables defined
nf_vars <- c(
    "projectDir",
    "fcid"
)
lapply(nf_vars, nf_var_check)

if (!exists("fcid")) {stop("'fcid' not defined!")}

# define data location
if (!exists("params.data_folder")) { # if data_loc not defined, use "data"
    data_loc="data"
} else {
    data_loc = params.data_folder
}

# TODO: Replace this code with non-functions.R code?

# run flow cell QC
step_seq_qc(fcid, data_loc)

# run index switching calculation
step_switching_calc(fcid, data_loc)

# stop(" *** stopped manually *** ") ##########################################