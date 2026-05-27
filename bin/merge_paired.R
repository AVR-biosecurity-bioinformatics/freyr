#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads_F                     <- args$reads_F
reads_R                     <- args$reads_R
seqs_F                      <- args$seqs_F
seqs_R                      <- args$seqs_R
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

### process variables 
if ( is.na(concat_unmerged) || concat_unmerged %in%  c("NA", "FALSE", "F")) {
    concat_unmerged <- FALSE
} else {
    concat_unmerged <- TRUE
}

### run R code
## process sample IDs to name the read and seq lists
## NOTE: the format of this list is different: "[a,b,c]" not "a b c"
sample_primers_list <-
    stringr::str_extract_all(
        sample_primers, 
        pattern = "[^\\s,\\[\\]]+" # extract all runs of characters that aren't spaces, commas or square brackets
        ) %>% 
    unlist()

## process F reads
reads_F_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_F, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()
names(reads_F_list) <- sample_primers_list

## process R reads
reads_R_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_R, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()
names(reads_R_list) <- sample_primers_list

## process F seqs
seqs_F_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_F, 
        pattern = "\\S+?_dada[1,2]F\\.rds" 
        ) %>% 
    unlist()

seqs_F_extracted <- lapply(seqs_F_list, readRDS)
names(seqs_F_extracted) <- sample_primers_list

## process R seqs
seqs_R_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_R, 
        pattern = "\\S+?_dada[1,2]R\\.rds" 
        ) %>% 
    unlist()

seqs_R_extracted <- lapply(seqs_R_list, readRDS)
names(seqs_R_extracted) <- sample_primers_list


# get names of NULL elements for F and R sequences
empty_seqs_F <- seqs_F_extracted[sapply(seqs_F_extracted, is.null)] %>% names
empty_seqs_R <- seqs_R_extracted[sapply(seqs_R_extracted, is.null)] %>% names

# remove empty sequence elements
seqs_F_pass <- seqs_F_extracted[!sapply(seqs_F_extracted, is.null)]
seqs_R_pass <- seqs_R_extracted[!sapply(seqs_R_extracted, is.null)]

# if either passing seq list is empty, throw error
if (length(seqs_F_pass) == 0 || length(seqs_R_pass) == 0){
    stop(paste0("\n***\nZero sequences made it through denoising for read group '",read_group,"' and primers '",primers,"'.\nConsider changing your read filtering parameters.\n***\n"))
}

# remove read files from list if there are NULL elements associated with them
reads_F_pass <- reads_F_list[!names(reads_F_list) %in% empty_seqs_F]
reads_R_pass <- reads_R_list[!names(reads_R_list) %in% empty_seqs_R]

# if either passing reads list is empty, throw error
if (length(reads_F_pass) == 0 || length(reads_R_pass) == 0){
    stop(paste0("\n***\nZero reads made it through denoising for read group '",read_group,"' and primers '",primers,"'.\nConsider changing your read filtering parameters.\n***\n"))
}

# vector of samples that pass (ie. have data)
sample_primers_pass <- sample_primers_list[!sample_primers_list %in% c(empty_seqs_F, empty_seqs_R)]
sample_primers_fail <- sample_primers_list[sample_primers_list %in% c(empty_seqs_F, empty_seqs_R)]

## merge pairs, keeping unmerged reads only if concat_unmerged is TRUE
mergers <- 
    dada2::mergePairs(
        dadaF = seqs_F_pass,
        derepF = reads_F_pass,
        dadaR = seqs_R_pass,
        derepR= reads_R_pass,
        verbose = TRUE,
        minOverlap = 12,
        trimOverhang = TRUE,
        returnRejects = concat_unmerged
    )


## reformat mergers as a list with one element if there is only one sample
if ( class(mergers) == "data.frame" ) {
    mergers <- list(mergers)
    names(mergers) <- sample_primers_pass
}

## TODO: write out unmerged reads? (pull from functions.R)

## concatenate unmerged reads
if ( concat_unmerged ) {
    message("concat_unmerged is set to TRUE - Concatenating unmerged forward and reverse reads")
    mergers_rescued <- mergers
    for(i in 1:length(mergers)) {
        if(any(!mergers[[i]]$accept)){
            # Get index of unmerged reads in table
            unmerged_index <- which(!mergers[[i]]$accept)
            # Get the forward and reverse reads for those reads
            unmerged_fwd <- seqs_F_pass[[i]]$sequence[mergers[[i]]$forward[unmerged_index]]
            unmerged_rev <- seqs_R_pass[[i]]$sequence[mergers[[i]]$reverse[unmerged_index]]
            
            unmerged_concat <- paste0(unmerged_fwd, "NNNNNNNNNN", rc(unmerged_rev))
            
            mergers_rescued[[i]]$sequence[unmerged_index] <- unmerged_concat
            mergers_rescued[[i]]$nmatch[unmerged_index] <- 0
            mergers_rescued[[i]]$nmismatch[unmerged_index] <- 0
            mergers_rescued[[i]]$nindel[unmerged_index] <- 0
            mergers_rescued[[i]]$prefer[unmerged_index] <- NA
            mergers_rescued[[i]]$accept[unmerged_index] <- TRUE
        } 
    }
    mergers <- mergers_rescued
}

# save number of merged reads per sample
getN <- function(x) sum(dada2::getUniques(x))

# create readsout output table
sapply(mergers, getN) %>% 
    as.data.frame() %>% 
    magrittr::set_colnames("pairs") %>%
    tibble::rownames_to_column(var = "sample_primers") %>%
    # add 0 abundance samples in with read counts of 0
    tibble::add_row(sample_primers = sample_primers_fail, pairs = 0) %>%
    dplyr::mutate(
        read_group = read_group,
        primers = primers,
        stage = "dada_mergereads"
    ) %>%
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>% 
    readr::write_csv(., paste0("dada_mergereads_",read_group,"_",primers,"_readsout.csv"))

# output mergers file as RDS
if (nrow(mergers[[1]]) > 0){
    saveRDS(mergers, "mergers.rds")
} else {
    saveRDS(NULL, "mergers.rds")
}

# stop(" *** stopped manually *** ") ##########################################

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})