#!/usr/bin/env Rscript

# check variables defined

### run R code
seqs_F <- readRDS(seqs_F)

seqs_R <- readRDS(seqs_R)

## merge pairs, keeping unmerged reads only if concat_unmerged is FALSE
if ( concat_unmerged ) {
    mergers <- dada2::mergePairs(
                        dadaF = seqs_F,
                        derepF = reads_F,
                        dadaR = seqs_R,
                        derepR= reads_R,
                        verbose = TRUE,
                        minOverlap = 12,
                        trimOverhang = TRUE,
                        returnRejects = TRUE
                    )
} else {
    mergers <- dada2::mergePairs(
                        dadaF = seqs_F,
                        derepF = reads_F,
                        dadaR = seqs_R,
                        derepR= reads_R,
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

saveRDS(mergers, paste0(sample_id,"_",pcr_primers,"_mergers.rds"))

# stop(" *** stopped manually *** ") ##########################################
