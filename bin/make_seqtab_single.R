#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads                       <- args$reads
seqs                        <- args$seqs
sample_primers              <- args$sample_primers
primers                     <- args$primers
read_group                  <- args$read_group
concat_unmerged             <- args$concat_unmerged

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

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
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))


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

# get names of NULL elements for sequences
empty_seqs <- seqs_extracted[sapply(seqs_extracted, is.null)] %>% names

# remove empty sequence elements
seqs_pass <- seqs_extracted[!sapply(seqs_extracted, is.null)]

# if either passing seq list is empty, throw error
if (length(seqs_pass) == 0 ){
    stop(paste0("\n***\nZero sequences made it through denoising for read group '",read_group,"' and primers '",primers,"'.\nConsider changing your read filtering parameters.\n***\n"))
}

# remove read files from list if there are NULL elements associated with them
reads_pass <- reads_list[!names(reads_list) %in% empty_seqs]

# if passing reads list is empty, throw error
if (length(reads_pass) == 0){
    stop(paste0("\n***\nZero reads made it through denoising for read group '",read_group,"' and primers '",primers,"'.\nConsider changing your read filtering parameters.\n***\n"))
}

# vector of samples that pass (ie. have data)
sample_primers_pass <- sample_primers_list[!sample_primers_list %in% empty_seqs]
sample_primers_fail <- sample_primers_list[sample_primers_list %in% empty_seqs]

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
    # add 0 abundance samples in with read counts of 0
    tibble::add_row(sample_primers = sample_primers_fail, pairs = 0) %>%
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

# add 0 abundance samples to seqtab
if (length(sample_primers_fail) > 0){
    seqtab_tibble[sample_primers_fail] <- 0L
}

# save DSS as .fasta
write_fasta(seq_DSS, file = paste0(read_group, "_", primers, "_seqs.fasta"))

# save seqtab_tidy as .csv
readr::write_csv(seqtab_tibble, paste0(read_group, "_", primers, "_seqtab_tibble.csv"))


# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})