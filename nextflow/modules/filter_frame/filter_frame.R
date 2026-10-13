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
coding                      <- args$coding
genetic_code                <- args$genetic_code

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

# Import phmm functions (formerly in taxreturn)
phmm_functions_file <- Sys.which("phmm_functions.R")
if (phmm_functions_file == "") {
    stop("Could not find functions.R on PATH") 
}

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

seqtab_list <- 
    seqtab_tibble_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) # read in seqtabs and store as list of tibbles

fasta_list <- 
    fasta_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

coding <-           as.logical(coding)

### run R code

# combine seqtabs into wide format
seqtab_combined <-
    seqtab_list %>%
    # pivot each tibble longer
    lapply(
        .,
        function(x){ # per tibble
            x %>%
            tidyr::pivot_longer(
                cols = !seq_name,
                names_to = "sample_primers",
                values_to = "abundance"
            )
        }
    ) %>%
    # bind tibbles together now columns all match
    dplyr::bind_rows() %>%
    # pivot wider, filling missing abundance with 0
    tidyr::pivot_wider(
        names_from = sample_primers,
        values_from = abundance, 
        values_fill = 0
    )

# combine sequences from .fasta files, removing redundancy, and convert to tibble (seq_name, sequence)
seqs_names <-  
    fasta_list %>%
    lapply(., as.character) %>%
    unlist(use.names = T) %>% 
    tibble::enframe(name = "seq_name", value = "sequence") %>%
    dplyr::group_by(seq_name, sequence) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

# check all sequences are unique
if( any(duplicated(seqs_names$sequence)) ){
    stop("Duplicated sequences found across input .fasta files!")
}

# check names in fastas are same as those in seqtabs
if ( !setequal(seqs_names$seq_name, seqtab_combined$seq_name) ){
    stop("Sequence names in seqtabs don't match sequence names in .fasta files!")
}

# create DSS object from tibble
seqs_dss <- 
    seqs_names %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

## frame checking
if(coding){
    # remove sequences with wrong frame  
    codon_filt_dss <- 
        codon_filter(
            seqs_dss, 
            genetic_code = genetic_code
        ) 
    
    # make tibble of name, sequence and whether it passed the PHMM filter
    out_tibble <- 
        seqs_names %>%
        dplyr::mutate(
            frame_filter = seq_name %in% names(codon_filt_dss)
        )
        
    } else {
    # all seqs pass filter
    out_tibble <- 
        seqs_names %>%
        dplyr::mutate(
            frame_filter = TRUE
        )
}

readr::write_csv(out_tibble, paste0(primers, "_frame_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})

#' Filter sequences containing stop codons
#'
#' @param x Sequences in DNAStringset or DNAbin format
#' @param genetic_code A genetic code for the Amino acid translation. set to 'SGC4' for Invertebrate mitochondrial or see all known codes at Biostrings::GENETIC_CODE_TABLE
#' @param tryrc Whether the reverse complemement should be evaluated if no frame without stop codons was found in the forward orientation.
#' @param resolve_draws How draws should be resolved when multiple possible frames produce sequences with no stop codons.
#' Options are "remove" to completely remove the sequence, or "majority" to pick the most common frame from the entire alignment.
#'
#' @return
#' @export
#' @importFrom ape as.DNAbin
#' @importFrom Biostrings reverseComplement
#' @importFrom DECIPHER RemoveGaps
#' @importFrom methods is
#'
codon_filter <- function(x, genetic_code = NULL, tryrc=TRUE, resolve_draws="majority"){
  if(is.null(genetic_code)){
    stop("genetic_code must not be NULL, set to 'SGC4' for Invertebrate mitochondrial or see Biostrings::GENETIC_CODE_TABLE for other options")
  }
  # Convert to DNAStringSet
  if (methods::is(x, "DNAbin")) {
    x <- DNAbin2DNAstringset(x, remove_gaps=FALSE)
    format <- "DNAbin"
  } else if(methods::is(x, "DNAStringSet")){
    format <- "DNAStringSet"
  } else {
    stop("x must be a DNAbin or DNAStringSet")
  }

  #Get reading frames
  frames <- get_reading_frame(DECIPHER::RemoveGaps(x), genetic_code = genetic_code, tryrc = tryrc, resolve_draws = resolve_draws)

  # Check if any sequences need RC (negative reading frame)
  to_rc <- sign(frames)==-1
  to_rc[is.na(to_rc)] <- FALSE
  if (any(to_rc)){
    message(sum(to_rc), " Sequences reverse complemented as no forward match was found")
    x[to_rc] <- Biostrings::reverseComplement(x[to_rc])
  }

  if(format=="DNAbin"){
    out <- ape::as.DNAbin(x[!is.na(frames)])
  } else if(format=="DNAStringSet"){
    out <- x[!is.na(frames)]
  }
  message(paste0(length(x) - length(out), " Sequences containing stop codons removed"))
  return(out)
}