#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

cpus                        <- as.integer(args$cpus)
primers                     <- args$primers
seqtab_tibble_list          <- args$seqtab_tibble_list
fasta_list                  <- args$fasta_list
minSampleFraction           <- as.numeric(args$minSampleFraction)


if (length(cpus) != 1L || is.na(cpus) || cpus < 1L) {
    stop("cpus must be a positive integer")
}

if (length(minSampleFraction) != 1L || !is.finite(minSampleFraction) || minSampleFraction < 0 || minSampleFraction > 1) {
    stop("minSampleFraction must be between 0 and 1")
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
    "dada2",
    "dplyr",
    "ggplot2",
    "patchwork",
    "readr",
    "stringr",
    "tibble",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

# Merge FASTA sequences first.
refseq <- fasta_list %>%
    stringr::str_split_1(pattern = " ") %>%
    lapply(Biostrings::readDNAStringSet) %>%
    merge_refseq()

# Read, merge, replace ASV IDs with DNA sequences, and transpose.
seqtab_matrix <- seqtab_tibble_list %>%
    stringr::str_split_1(pattern = " ") %>%
    lapply(function(file) {
        readr::read_csv(
            file,
            col_types = readr::cols(
                seq_name = readr::col_character(),
                .default = readr::col_integer()
            ),
            show_col_types = FALSE
        ) %>%
            tibble::column_to_rownames(var = "seq_name") %>%
            as.matrix()
    }) %>%
    merge_seqtab(taxa_are_rows = FALSE, refseq = refseq)

# remove chimeras
seqtab_nochim <- 
    dada2::removeBimeraDenovo(
  		seqtab_matrix, 
  		method = "consensus",
        minSampleFraction = minSampleFraction,
        multithread = cpus
    )

# Output sequence IDs, DNA sequences, and filter results.
out_tibble <- tibble::tibble(
    seq_name = names(refseq),
    sequence = as.character(refseq),
    chimera_filter = sequence %in% colnames(seqtab_nochim)
)

readr::write_csv(out_tibble,paste0(primers, "_chimera_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})