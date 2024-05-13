#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    # "Biostrings",
    # "bs4Dash",
    # "clustermq",
    # "dada2",
    # "DECIPHER",
    "dplyr",
    # "future",
    "ggplot2",
    # "gridExtra",
    # "gt",
    "magrittr",
    # "markdown",
    # "ngsReports",
    # "patchwork",
    # "phyloseq",
    # "pingr",
    "purrr",
    "readr",
    # "rlang",
    # "rstudioapi",
    "savR",
    # "scales",
    "seqateurs",
    # "shiny",
    # "shinybusy",
    # "shinyWidgets",
    # "ShortRead",
    "stringr",
    # "taxreturn",
    "tibble",
    "tidyr",
    # "vegan",
    # "visNetwork",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE))

#### TODO: Pull stats from top of index_switching.pdf and print to a final run report
### also flag index combos that have higher than average + (some stat)

#### TODO: use index_switch_calc.txt in jack_notes to run same process but in bash, which should be much faster

# check variables defined
if (!exists("fcid")) {stop("'fcid' not defined!")}

# define data location
if (!exists("params.data_folder")) { # if data_loc not defined, use "data"
    data_loc="data"
} else {
    data_loc = params.data_folder
}

# TODO: Replace this code with non-functions.R code?

# run flow cell QC
step_seq_qc(fcid)

# run index switching calculation
step_switching_calc(fcid)

stop(" *** stopped manually *** ") ##########################################