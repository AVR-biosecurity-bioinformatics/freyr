#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dada2",
    "dplyr",
    "magrittr",
    "readr",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "read_group",
    "primers",
    "sample_primers",
    "reads",
    "seqs",
    "concat_unmerged"
)
lapply(nf_vars, nf_var_check)

### run R code
## process sample IDs to name the read and seq lists
## NOTE: the format of this list is different: "[a,b,c]" not "a b c"
sample_primers_list <-
    stringr::str_extract_all(
        sample_primers, 
        pattern = "[^\\s,\\[\\]]+" # extract all runs of characters that aren't spaces, commas or square brackets
        ) %>% 
    unlist()

## process single-end reads
reads_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()
names(reads_list) <- sample_primers_list

## process single-end seqs
seqs_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs, 
        pattern = "\\S+?_dada[1,2]S\\.rds" 
        ) %>% 
    unlist()

seqs_extracted <- lapply(seqs_list, readRDS)
names(seqs_extracted) <- sample_primers_list

## reformat seqs_extracted as a list with one element if there is only one sample
if ( class(seqs_extracted) == "data.frame" ) {
  seqs_extracted <- list(seqs_extracted)
  names(seqs_extracted) <- sample_primers_list
}

## make sequence table
seqtab <- dada2::makeSequenceTable(seqs_extracted)

saveRDS(seqtab, paste0(read_group,"_",primers,"_seqtab.rds"))

### TODO: recode this
# save number of merged reads per sample
getN <- function(x) sum(dada2::getUniques(x))

sapply(seqs_extracted, getN) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames("pairs") %>%
    tibble::rownames_to_column(var = "sample_primers") %>%
    dplyr::mutate(
        read_group = read_group,
        primers = primers,
        stage = "make_seqtab"
    ) %>%
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>% 
    readr::write_csv(., paste0("dada_mergereads_",read_group,"_",primers,"_readsout.csv"))


##### make new outputs ----------------------------------------------------------------------------


### convert DADA2 format sequence table to both .fasta of sequences with hashes and .csv of sample abundances

# sequences as vector
seq_vec <- seqtab %>% dada2::getSequences()

# sequences as hash
hash_vec <- seq_vec %>% lapply(., rlang::hash) %>% unlist()

# named vector of sequences
names(seq_vec) <- hash_vec

# named DSS
seq_DSS <- Biostrings::DNAStringSet(seq_vec)

# tibble of hash and seq
seq_tibble <- tibble::enframe(seq_vec, name = "seq_name", value = "sequence")

# convert seqtab to one row per sequence, ID with hash
seqtab_tibble <- 
    seqtab %>% 
    tidyr::as_tibble(rownames = "sample_primers") %>%
    tidyr::pivot_longer(cols = !sample_primers, names_to = "sequence", values_to = "abundance") %>%
    # join to seq_tibble
    dplyr::left_join(., seq_tibble, by = "sequence") %>%
    dplyr::select(-sequence) %>% # remove sequence
    tidyr::pivot_wider(names_from = sample_primers, values_from = abundance)


# save DSS as .fasta
write_fasta(seq_DSS, file = paste0(read_group, "_", primers, "_seqs.fasta"))

# save seqtab_tidy as .csv
readr::write_csv(seqtab_tibble, paste0(read_group, "_", primers, "_seqtab_tibble.csv"))


# stop(" *** stopped manually *** ") ##########################################
