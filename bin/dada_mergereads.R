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
    "reads_F",
    "reads_R",
    "seqs_F",
    "seqs_R",
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

## process F reads
reads_F_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_F, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()
names(reads_F_list) <- sample_id_list

## process R reads
reads_R_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_R, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()
names(reads_R_list) <- sample_id_list

## process F seqs
seqs_F_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_F, 
        pattern = "\\S+?_dada[1,2]F\\.rds" 
        ) %>% 
    unlist()

seqs_F_extracted <- lapply(seqs_F_list, readRDS)
names(seqs_F_extracted) <- sample_id_list

## process R seqs
seqs_R_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_R, 
        pattern = "\\S+?_dada[1,2]R\\.rds" 
        ) %>% 
    unlist()

seqs_R_extracted <- lapply(seqs_R_list, readRDS)
names(seqs_R_extracted) <- sample_id_list

## merge pairs, keeping unmerged reads only if concat_unmerged is FALSE
if ( concat_unmerged ) {
    mergers <- dada2::mergePairs(
        dadaF = seqs_F_extracted,
        derepF = reads_F_list,
        dadaR = seqs_R_extracted,
        derepR= reads_R_list,
        verbose = TRUE,
        minOverlap = 12,
        trimOverhang = TRUE,
        returnRejects = TRUE
    )
} else {
    mergers <- dada2::mergePairs(
        dadaF = seqs_F_extracted,
        derepF = reads_F_list,
        dadaR = seqs_R_extracted,
        derepR= reads_R_list,
        verbose = TRUE,
        minOverlap = 12,
        trimOverhang = TRUE,
        returnRejects = FALSE
    )
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
            unmerged_fwd <- seqs_F[[i]]$sequence[mergers[[i]]$forward[unmerged_index]]
            unmerged_rev <- seqs_R[[i]]$sequence[mergers[[i]]$reverse[unmerged_index]]
            
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

# Construct sequence table for fcid x pcr_primers from merged reads per sample
seqtab <- dada2::makeSequenceTable(mergers)

saveRDS(seqtab, paste0(fcid,"_",pcr_primers,"_seqtab.rds"))

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
