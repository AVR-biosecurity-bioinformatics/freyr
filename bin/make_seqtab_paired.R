#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

mergers_list                <- args$mergers_list
sample_primers              <- args$sample_primers
primers                     <- args$primers
read_group                  <- args$read_group

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

### process variables 
mergers_list <- 
  mergers_list %>% 
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
  ) 

### run R code
# get names of NULL elements for F and R sequences
empty_mergers <- mergers_list[sapply(mergers_list, is.null)] %>% names

# remove empty merger elements
mergers_pass <- mergers_list[!sapply(mergers_list, is.null)]

if (length(mergers_pass) == 0){
  stop("All samples have no sequences")
}

seqtab <- dada2::makeSequenceTable(mergers_pass)

seq_vec <- seqtab %>% dada2::getSequences()

# sequences as hash
hash_vec <- seq_vec %>% lapply(., digest::digest, algo = "sha256") %>% unlist()

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
if (length(empty_mergers) > 0){
    seqtab_tibble[empty_mergers] <- 0L
}

##### make new outputs ----------------------------------------------------------------------------


### convert DADA2 format sequence table to both .fasta of sequences with hashes and .csv of sample abundances

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