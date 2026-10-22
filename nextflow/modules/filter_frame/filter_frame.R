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
coding                      <- as.logical(args$coding)
genetic_code                <- args$genetic_code

if (length(coding) != 1L || is.na(coding)) {
    stop("coding must be TRUE or FALSE")
}

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

# Import phmm functions (formerly in taxreturn)
phmm_functions_file <- Sys.which("phmm_functions.R")
if (phmm_functions_file == "") {
    stop("Could not find phmm_functions.R on PATH") 
}
sys.source(phmm_functions_file, envir = .GlobalEnv)

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

# Merge FASTA sequences
seqs_dss <- fasta_list %>%
    stringr::str_split_1(pattern = " ") %>%
    lapply(Biostrings::readDNAStringSet) %>%
    merge_refseq()

# Default to passing all sequences for non-coding markers.
frame_filter <- rep(TRUE, length(seqs_dss))

## frame checking - only with coding
if(coding){
    # remove sequences with wrong frame  
    codon_filt_dss <- 
        codon_filter(
            seqs_dss, 
            genetic_code = genetic_code
        ) 
    
    frame_filter <- names(seqs_dss) %in% names(codon_filt_dss)      
 } 

# TRUE means the sequence passed frame filtering.
out_tibble <- tibble::tibble(
    seq_name = names(seqs_dss),
    sequence = as.character(seqs_dss),
    frame_filter = frame_filter
)

readr::write_csv(out_tibble, paste0(primers, "_frame_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})
