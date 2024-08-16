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
    "fcid",
    "pcr_primers",
    "sample_id",
    "reads",
    "seqs",
    "concat_unmerged"
)
lapply(nf_vars, nf_var_check)

### run R code
## process sample IDs to name the read and seq lists
## NOTE: the format of this list is different: "[a,b,c]" not "a b c"
sample_id_list <-
    stringr::str_extract_all(
        sample_id, 
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
names(reads_list) <- sample_id_list

## process single-end seqs
seqs_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs, 
        pattern = "\\S+?_dada[1,2]F\\.rds" 
        ) %>% 
    unlist()

seqs_extracted <- lapply(seqs_list, readRDS)
names(seqs_extracted) <- sample_id_list

## reformat seqs_extracted as a list with one element if there is only one sample
if ( class(seqs_extracted) == "data.frame" ) {
  seqs_extracted <- list(seqs_extracted)
  names(seqs_extracted) <- sample_id_list
}

## make sequence table
seqtab <- dada2::makeSequenceTable(seqs_extracted)

saveRDS(seqtab, paste0(fcid,"_",pcr_primers,"_seqtab.rds"))

### TODO: recode this
# save number of merged reads per sample
getN <- function(x) sum(dada2::getUniques(x))

sapply(mergers, getN) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames("pairs") %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    dplyr::mutate(
        fcid = fcid,
        pcr_primers = pcr_primers,
        stage = "dada_mergereads"
    ) %>%
    dplyr::select(stage, sample_id, fcid, pcr_primers, pairs) %>% 
    readr::write_csv(., paste0("dada_mergereads_",fcid,"_",pcr_primers,"_readsout.csv"))


# stop(" *** stopped manually *** ") ##########################################
