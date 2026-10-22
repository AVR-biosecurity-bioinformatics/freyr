#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
seqtab_tibble_list          <- args$seqtab_tibble_list
fasta_list                  <- args$fasta_list
asv_min_length              <- suppressWarnings(as.numeric(args$asv_min_length))
asv_max_length              <- suppressWarnings(as.numeric(args$asv_max_length))

 # NA disables the corresponding bound.
if (is.na(asv_min_length)) {
    asv_min_length <- -Inf
}
if (is.na(asv_max_length)) {
    asv_max_length <- Inf
}

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "readr",
    "stringr",
    "tibble",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

## check and define variables

# Merge FASTA sequences
seqs_dss <- fasta_list %>%
    stringr::str_split_1(pattern = " ") %>%
    lapply(Biostrings::readDNAStringSet) %>%
    merge_refseq()

# Calculate lengths without converting sequences to character strings.
seq_lengths <- BiocGenerics::width(seqs_dss)

# TRUE means the sequence passed length filtering.
out_tibble <- tibble::tibble(
    seq_name = names(seqs_dss),
    sequence = as.character(seqs_dss),
    length_filter = seq_lengths > asv_min_length &
                    seq_lengths < asv_max_length
)

readr::write_csv(out_tibble, paste0(primers, "_length_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})