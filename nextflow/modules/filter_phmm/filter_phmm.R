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
coding                      <- args$coding

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

seqs <- Biostrings::readDNAStringSet(fasta)

coding <-           as.logical(coding)

if(is.na(phmm))             {phmm <- NULL}
if(is.na(for_primer_seq))   {for_primer_seq <- NULL}
if(is.na(rev_primer_seq))   {rev_primer_seq <- NULL}

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


primer_seqs <- c(for_primer_seq, rev_primer_seq)

### run R code

# combine sequences from .fasta files, removing redundancy, and convert to tibble (seq_name, sequence)
seqs_names <-  
    seqs %>%
    as.character() %>%
    tibble::enframe(name = "seq_name", value = "sequence") %>%
    dplyr::group_by(seq_name, sequence) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

# check all sequences are unique
if( any(duplicated(seqs_names$sequence)) ){
    stop("Duplicated sequences found across input .fasta files!")
}

# create DSS object from tibble
seqs_dss <- 
    seqs_names %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

## subset PHMM if primers were provided
if (is(phmm_model, "PHMM") && !is.null(primer_seqs)){
    # Check that one of the two primers can bind
    Fbind <- get_binding_position(primer_seqs[1], model = phmm_model, tryRC = TRUE, min_score = 10)
    Rbind <- get_binding_position(primer_seqs[2], model = phmm_model, tryRC = TRUE, min_score = 10)
    if(!is.na(Fbind$start) & !is.na(Rbind$start)){
        phmm_model <- subset_model(phmm_model, primers = primer_seqs)
    } else if(!is.na(Fbind$start) & is.na(Rbind$start)){
        # Reverse primer not found - Try with subsets
        for(r in seq(1, nchar(primer_seqs[2])-10, 1)){ #Minimum length of 10 as this has to match minscore
            Rbind <- get_binding_position(
                    str_remove(primer_seqs[2], paste0("^.{1,",r,"}")), 
                    model = phmm_model, 
                    tryRC = TRUE, 
                    min_score = 10
                )
            if (!is.na(Rbind$start)) {
                primer_seqs[2] <- stringr::str_remove(primer_seqs[2], paste0("^.{1,",r,"}"))
                break
            }
        }
        phmm_model <- subset_model(phmm_model, primers = primer_seqs)
    } else  if(is.na(Fbind$start) & !is.na(Rbind$start)){
        # Forward primer not found - Try with subsets
        for(r in seq(1, nchar(primer_seqs[1])-10, 1)){ #Minimum length of 10 as this has to match minscore
            Rbind <- get_binding_position(
                    stringr::str_remove(primer_seqs[1], paste0("^.{1,",r,"}")), 
                    model = phmm_model, 
                    tryRC = TRUE, 
                    min_score = 10
                )
            if (!is.na(Rbind$start)) {
                primer_seqs[1] <- stringr::str_remove(primers[1], paste0("^.{1,",r,"}"))
                break
            }
        }
        phmm_model <- subset_model(phmm_model, primers = primer_seqs)
    }
}

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
    
    # make tibble of name, sequence and whether it passed the PHMM filter
    out_tibble <- 
        seqs_names %>%
        dplyr::mutate(
            phmm_filter = seq_name %in% names(phmm_filt_seqtab)
        )

} else {
    # all seqs pass filter
    out_tibble <- 
        seqs_names %>%
        dplyr::mutate(
            phmm_filter = TRUE
        )
}

readr::write_csv(out_tibble, paste0(primers, "_phmm_filter.csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})

#' Subset PHMM
#' Modified from Cameron Nugent's COIL package https://github.com/CNuge/coil/blob/master/R/subset_PHMM.r
#' @param x A profile hidden Markov model (a "PHMM" object) generated by the aphid R package.
#' @param primers A character vector of length 2 contaning the forward and reverse primer sequences.
#' @param positions A numeric vector of length 2 containing the start and end position for model trimming. Used only if primers = NULL.
#' @param trimprimers logical indicating whether the primer-binding sites should be removed from the sequences in the returned model. Used only if positions = NULL.
#' @param min_score Minimum score for the viterbi alignment.
#' @param quiet Whether to print progress to console.
#'
#' @return
#' @export
#'
#' @examples
subset_model <- function(x, primers=NULL, positions=NULL, trimprimers=FALSE, min_score=10, quiet=FALSE){

  # Check inputs
  if (is.numeric(positions) & length(positions) == 2){
    start_bind <- sort(positions)[1]
    end_bind <- sort(positions)[2]
  } else if (is.null(positions) & (is.character(primers) & length(primers)== 2)){
    Fprimer <- primers[1]
    Rprimer <- primers[2]

    Fbind <-  get_binding_position(Fprimer, model=x, tryRC=TRUE, quiet=quiet, min_score = min_score)
    Rbind <-  get_binding_position(Rprimer, model=x, tryRC=TRUE, quiet=quiet, min_score = min_score)
    if(isTRUE(trimprimers)){
      start_bind <- Fbind$end + 1
      end_bind <- Rbind$start - 1
      if(!quiet){message("Primer binding positions will be trimmed from model")}
    } else if(isFALSE(trimprimers)){
      start_bind <- Fbind$start
      end_bind <- Rbind$end
    } else{
      stop("trimprimers must be either TRUE or FALSE")
    }

  } else if (is.null(positions) & is.null(primers)){
    stop("A value must be provided to either positions or primers")
  } else if(is.null(primers) & (!is.numeric(positions) | !length(positions) == 2) ){
    stop("positions must be a numeric vector of length 2 containing the start and end position for model trimming")
  } else if(is.null(positions) & (!is.character(primers) | !length(primers) == 2) ){
    stop("primers must be a character vector of length 2, with the forward primer as the first element, and reverse primer as the second element")
  } else if(!is.null(primers) & !is.null(positions)){
    stop("Only one of primers or positions can be entered at a time")
  } else {
    stop("A value must be provided to either positions or primers")
  }
  if(is.na(start_bind) & is.numeric(end_bind)){
    stop(paste0("Could not find binding position for Forward primer: ", Fprimer))
  } else if(is.na(end_bind) & is.numeric(start_bind)){
    stop(paste0("Could not find binding position for Reverse primer: ", Rprimer))
  } else if(is.na(start_bind) & is.na(end_bind)){
    stop(paste0("Could not find binding position for Forward primer: ", Fprimer, " and Reverse primer: ", Rprimer))
  }
  # Check that positions are right
  if(end_bind < start_bind){
    stop("Index error, end position must match or exceed the start position.")
  }
  if(end_bind > x$size){
    stop(paste0("Index error, end position out of bounds. Input PHMM only has a length of: ", x$size))
  }

  #get the map positions for subsetting the alignment-based fields
  map_start <- x$map[[start_bind]]
  map_end <- x$map[[end_bind]]

  #subset to get the new inserts
  new_inserts <- x$inserts[map_start:map_end]

  #subset to get the new insert lengths
  new_insert_lengths <- x$insertlengths[map_start:map_end]
  #reindex the insert length labels
  if(!is.null(new_insert_lengths)){
    names(new_insert_lengths) <- 0:(length(new_insert_lengths)-1)
  }
  #subset the mask
  new_mask <- x$mask[map_start:map_end]
  #subset the map
  new_map <- x$map[start_bind:end_bind]

  #set the first remaining map entry equal to one,
  #make the necessary relative adjustment to all other positions
  new_map <- new_map-(new_map[1]-1)
  #reindex the labels
  if(!is.null(new_map)){
    names(new_map) <- 1:(length(new_map))
  }

  #Subset the A and E matrixies, reindex the positions
  #A is 0 indexed, E is 1 indexed
  new_A <- x$A[,start_bind:end_bind]
  colnames(new_A) <- 0:(length(colnames(new_A))-1)
  new_E <- x$E[,start_bind:end_bind]
  colnames(new_E) <- 1:length(colnames(new_A))

  #Create the new PHMM object
  newPHMM = structure(list(name = x$name,
                           description = x$description,
                           size = (end_bind - start_bind)+1,
                           alphabet = x$alphabet,
                           A = new_A,
                           E = new_E,
                           qa = x$qa,
                           qe = x$qe,
                           inserts = new_inserts,
                           insertlengths = new_insert_lengths,
                           map = new_map,
                           date = x$date,
                           nseq = x$nseq,
                           weights = x$weights,
                           reference = x$reference,
                           mask = new_mask
  ), class =  "PHMM")
  message(paste0("PHMM model subset to ", newPHMM$size, " base pairs"))
  return(newPHMM)
}