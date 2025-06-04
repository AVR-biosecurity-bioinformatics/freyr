#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    # "dada2",
    "dplyr",
    # "ggplot2",
    # "patchwork",
    "readr",
    "stringr",
    "taxreturn",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "pcr_primers",
    "seqtab_tibble_list",
    "fasta_list",
    "phmm",
    "for_primer_seq",
    "rev_primer_seq",
    "coding"
)
lapply(nf_vars, nf_var_check)

## check and define variables

seqtab_list <- 
    seqtab_tibble_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) # read in seqtabs and store as list of tibbles

fasta_list <- 
    fasta_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

phmm <-             parse_nf_var_repeat(phmm)
for_primer_seq <-   parse_nf_var_repeat(for_primer_seq)
rev_primer_seq <-   parse_nf_var_repeat(rev_primer_seq)
coding <-           parse_nf_var_repeat(coding) %>% as.logical

if(is.na(phmm))             {phmm <- NULL}
if(is.na(for_primer_seq))   {for_primer_seq <- NULL}
if(is.na(rev_primer_seq))   {rev_primer_seq <- NULL}

primers <- c(for_primer_seq, rev_primer_seq)

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
                names_to = "sample_id",
                values_to = "abundance"
            )
        }
    ) %>%
    # bind tibbles together now columns all match
    dplyr::bind_rows() %>%
    # pivot wider, filling missing abundance with 0
    tidyr::pivot_wider(
        names_from = sample_id,
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

## Load in profile hidden markov model if provided
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

## subset PHMM if primers were provided
if (is(phmm_model, "PHMM") && !is.null(primers)){
    # Check that one of the two primers can bind
    Fbind <- taxreturn::get_binding_position(primers[1], model = phmm_model, tryRC = TRUE, min_score = 10)
    Rbind <- taxreturn::get_binding_position(primers[2], model = phmm_model, tryRC = TRUE, min_score = 10)
    if(!is.na(Fbind$start) & !is.na(Rbind$start)){
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else if(!is.na(Fbind$start) & is.na(Rbind$start)){
        # Reverse primer not found - Try with subsets
        for(r in seq(1, nchar(primers[2])-10, 1)){ #Minimum length of 10 as this has to match minscore
            Rbind <- get_binding_position(
                str_remove(primers[2], paste0("^.{1,",r,"}")), 
                model = phmm_model, 
                tryRC = TRUE, 
                min_score = 10
                )
            if (!is.na(Rbind$start)) {
                primers[2] <- stringr::str_remove(primers[2], paste0("^.{1,",r,"}"))
                break
            }
        }
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else  if(is.na(Fbind$start) & !is.na(Rbind$start)){
        # Forward primer not found - Try with subsets
        for(r in seq(1, nchar(primers[1])-10, 1)){ #Minimum length of 10 as this has to match minscore
        Rbind <- taxreturn::get_binding_position(stringr::str_remove(primers[1], paste0("^.{1,",r,"}")), model = phmm_model, tryRC = TRUE, min_score = 10)
        if (!is.na(Rbind$start)) {
            primers[1] <- stringr::str_remove(primers[1], paste0("^.{1,",r,"}"))
            break
        }
        }
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    }
}

## PHMM filtering
if (is(phmm_model, "PHMM")){
    # remove sequences that don't align to PHMM
    phmm_filt_seqtab <- 
        taxreturn::map_to_model(
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

readr::write_csv(out_tibble, paste0(pcr_primers, "_phmm_filter.csv"))
