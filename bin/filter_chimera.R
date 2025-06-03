#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "patchwork",
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
    "minSampleFraction"
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

minSampleFraction <- as.numeric(minSampleFraction)

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

# make a dada2-style seqtab from both tibbles
seqtab_matrix <- 
    seqtab_combined %>%
    # add seq string
    dplyr::left_join(., seqs_names, by = "seq_name") %>%
    # remove seq_name
    dplyr::select(-seq_name) %>%
    # pivot longer
    tidyr::pivot_longer(cols = !sequence, names_to = "sample_id", values_to = "abundance") %>%
    dplyr::mutate(abundance = as.integer(abundance)) %>%
    # pivot wider 
    tidyr::pivot_wider(names_from = sequence, values_from = abundance) %>%
    # sample_id as rownames
    tibble::column_to_rownames(var = "sample_id") %>%
    # convert to matrix
    as.matrix() 

# remove chimeras
seqtab_nochim <- 
    dada2::removeBimeraDenovo(
  		seqtab_matrix, 
  		method = "consensus",
        minSampleFraction = minSampleFraction
    )

# output table of which sequences passed or failed filter
out_tibble <- 
  seqs_names %>%
	dplyr::mutate(chimera_filter = sequence %in% colnames(seqtab_nochim))

readr::write_csv(out_tibble, paste0(pcr_primers, "_chimera_filter.csv"))

# output read tracking tibble (to be modified in the "merge" process to add fcid column)
read_tracking_out <- 
    out_tibble %>%
    dplyr::left_join(., seqtab_combined, by = "seq_name") %>%
    dplyr::rename("filter_chimera" = chimera_filter) %>%
    tidyr::pivot_longer(
        cols = !c(seq_name, sequence, filter_chimera),
        names_to = "sample_id",
        values_to = "abundance"
    ) %>%
    dplyr::select(-seq_name, -sequence) %>%
    dplyr::mutate(pcr_primers = pcr_primers) %>%
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("filter_"), 
        names_to = "stage",
        values_to = "pass"
    ) %>%
    dplyr::group_by(sample_id, pcr_primers, stage, pass) %>%
    dplyr::summarise(pairs = sum(abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(pass == TRUE) %>%
    dplyr::select(-pass) %>%
    dplyr::select(stage, sample_id, pcr_primers, pairs)

readr::write_csv(read_tracking_out, paste0("filter_chimeras_",pcr_primers,"_readsout.csv"))

