#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

mergers_paths               <- args$mergers_paths
sample_primers              <- args$sample_primers
primers                     <- args$primers
read_group                  <- args$read_group

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dada2",
    "dplyr",
    "magrittr",
    "readr",
    "stringr",
    "tibble",
    "digest",
    "Biostrings",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

### process variables 
mergers_list <- 
  mergers_paths %>% 
  # extract all runs of characters that aren't spaces, commas or square brackets
  stringr::str_extract_all(
    ., 
    pattern = "[^\\s,\\[\\]]+"
    ) %>% 
  unlist() %>% 
  sapply(
    ., 
    function(x){
      readRDS(x)
    },
    USE.NAMES = F
  ) #%>%
  #unlist(., recursive = F)

### run R code
# get names of NULL elements for F and R sequences
empty_mergers <- mergers_list[sapply(mergers_list, is.null)] %>% names

# remove empty merger elements
mergers_pass <- mergers_list[!sapply(mergers_list, is.null)]

if (length(mergers_pass) == 0){
  stop("All samples have no sequences")
}

### convert DADA2 format sequence table to both .fasta of sequences with hashes and .csv of sample abundances

seqtab <- dada2::makeSequenceTable(mergers_pass)

seq_vec <- dada2::getSequences(seqtab)

# sequences as hash
hash_vec <- vapply(
    seq_vec,
    digest::digest,
    character(1),
    algo = "sha256"
)

# Replace sequence columns with hashes
colnames(seqtab) <- hash_vec

# Export one row per sequence: seq_name followed by sample abundances.
# Do not include a sequence column: downstream expects all other
# columns to contain numeric sample counts.
seqtab_out <- data.frame(
    seq_name = hash_vec,
    t(seqtab),
    check.names = FALSE,
    row.names = NULL
)

# Add samples with no sequences
if (length(empty_mergers) > 0) {
    seqtab_out[empty_mergers] <- 0L
}

# Write outputs
seq_DSS <- Biostrings::DNAStringSet(seq_vec)
names(seq_DSS) <- hash_vec
write_fasta(seq_DSS, file = paste0(read_group, "_", primers, "_seqs.fasta"))

data.table::fwrite(seqtab_out, paste0(read_group, "_", primers, "_seqtab_tibble.csv"))

# stop(" *** stopped manually *** ") ##########################################

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})