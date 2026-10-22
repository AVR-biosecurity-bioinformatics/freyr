#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
fasta                       <- args$fasta
phmm                        <- args$phmm
for_primer_seq              <- args$for_primer_seq
rev_primer_seq              <- args$rev_primer_seq
coding                      <- as.logical(args$coding)

if(is.na(phmm))             {phmm <- NULL}
if(is.na(for_primer_seq))   {for_primer_seq <- NULL}
if(is.na(rev_primer_seq))   {rev_primer_seq <- NULL}

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
    "aphid",
    "future",
    "furrr",
    "purrr",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

## check and define variables

# Read sequences directly; assume unique, validated ASVs upstream.
seqs_dss <- Biostrings::readDNAStringSet(fasta)

# Load in profile hidden markov model if provided
if(is.character(phmm) && stringr::str_detect(phmm, ".rds")){
    phmm_model <- readRDS(phmm)
    message("Loaded PHMM from file")
} else if (is(phmm, "PHMM")){
    phmm_model <- phmm
    message("Loaded PHMM from R object")
} else {
    phmm_model <- NULL
    message("Running analysis with no PHMM")
}


# Subset PHMM model only when both primers are available and can be located.
if (methods::is(phmm_model, "PHMM") && !is.null(for_primer_seq) && !is.null(rev_primer_seq)) {
    forward <- find_binding_primer(for_primer_seq, phmm_model)
    reverse <- find_binding_primer(rev_primer_seq, phmm_model)
    
    # Preserve the original strategy: trim a missing primer only
    # when the other full-length primer binds.
    if (!is.null(forward) && is.null(reverse)) {
        reverse <- find_binding_primer(
            rev_primer_seq,
            phmm_model,
            try_trim = TRUE
        )
    } else if (is.null(forward) && !is.null(reverse)) {
        forward <- find_binding_primer(
            for_primer_seq,
            phmm_model,
            try_trim = TRUE
        )
    }
    if (!is.null(forward) && !is.null(reverse)) {
        phmm_model <- subset_model(
            phmm_model,
            primers = c(forward, reverse)
        )
    } else {
        warning(
            "Could not locate both primers; using the full PHMM."
        )
    }
}

# Default to passing all sequences when no PHMM is supplied.
phmm_filter <- rep(TRUE, length(seqs_dss))

## PHMM filtering
if (is(phmm_model, "PHMM")){
    # remove sequences that don't align to PHMM
    set.seed(1); phmm_filt_seqtab <- 
        map_to_model(
            seqs_dss, 
            model = phmm_model, 
            min_score = 100, 
            min_length = 100,
            shave = FALSE, 
            check_frame = coding, 
            kmer_threshold = 0.5, 
            k = 5, 
            extra = "fill"
        )
    
      phmm_filter <- names(seqs_dss) %in% names(phmm_filt_seqtab)

}

# TRUE means the sequence passed PHMM filtering.
out_tibble <- tibble::tibble(
    seq_name = names(seqs_dss),
    sequence = as.character(seqs_dss),
    phmm_filter = phmm_filter
)
readr::write_csv(out_tibble, paste0(primers, "_phmm_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})
