#!/usr/bin/env Rscript

# check variables defined
if ( !concat_unmerged %in% c(TRUE, FALSE) ) {
    stop("**** 'concat_unmerged' must be TRUE or FALSE ****")
}

### run R code
## process F reads
reads_F_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_F, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()

## process R reads
reads_R_list <- # convert input reads list from Groovy format to R format
    stringr::str_extract_all(
        reads_R, 
        pattern = "\\S+?\\.fastq\\.gz|\\S+?\\.fastq|\\S+?\\.fq\\.gz|\\S+?\\.fq" 
        ) %>% 
    unlist()

## process F seqs
seqs_F_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_F, 
        pattern = "\\S+?_dada[1,2]F\\.rds" 
        ) %>% 
    unlist()

seqs_F_extracted <- list() # new list
for (i in 1:length(seqs_F_list)) { # loop through reading .rds files and add to list
    seq <- readRDS(seqs_F_list[i])
    seqs_F_extracted <- append(seqs_F_extracted, seq)
}

print(seqs_F_extracted[1])
print(seqs_F_extracted[2])

stop(" *** stopped manually *** ") ##########################################


## process R seqs
seqs_R_list <- # convert input sequences list from Groovy format to R format
    stringr::str_extract_all(
        seqs_R, 
        pattern = "\\S+?_dada[1,2]R\\.rds" 
        ) %>% 
    unlist()

seqs_R_extracted <- list() # new list
for (i in 1:length(seqs_R_list)) { # loop through reading .rds files and add to list
    seq <- readRDS(seqs_R_list[i])
    seqs_R_extracted <- append(seqs_R_extracted, seq)
}





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

print(mergers)



# Construct sequence table for fcid x pcr_primers from merged reads per sample
seqtab <- dada2::makeSequenceTable(mergers)

print(seqtab)



saveRDS(seqtab, paste0(sample_id,"_",pcr_primers,"_seqtab.rds"))

# # Track reads
# getN <- function(x) sum(getUniques(x))
# res <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) %>%
#     magrittr::set_colnames(c("dadaFs", "dadaRs", "merged")) %>%
#     as.data.frame() %>%
#     rownames_to_column("sample_id") %>%
#     dplyr::mutate(sample_id = stringr::str_remove(basename(sample_id), pattern="_S[0-9]+_R[1-2]_.*$")) %>%
#     as_tibble()
# return(res)

