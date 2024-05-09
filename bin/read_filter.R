#!/usr/bin/env Rscript

# check variables defined


### run R code

## using dada2::filterAndTrim directly
res <- dada2::filterAndTrim(
    fwd = fwd_reads, 
    filt = paste0(sample_id,"_",target_gene,"_",pcr_primers,"_filter_R1.fastq.gz"), # gets saved to working dir
    rev = rev_reads, 
    filt.rev = paste0(sample_id,"_",target_gene,"_",pcr_primers,"_filter_R2.fastq.gz"), # gets saved to working dir
    minLen = as.numeric(read_min_length), 
    maxLen = as.numeric(read_max_length), 
    maxEE = as.numeric(read_max_ee), 
    truncLen = as.numeric(read_trunc_length),
    trimLeft = as.numeric(read_trim_left), 
    trimRight = as.numeric(read_trim_right), 
    rm.phix = TRUE, 
    multithread = FALSE, 
    compress = TRUE, 
    verbose = FALSE
)

## extract output read counts and save to file
reads_out <- res %>%
    tibble::as_tibble() %>%
    pull(reads.out)

c("read_filter", sample_id, fcid, pcr_primers, reads_out, reads_out) %>% 
    tibble::as_tibble() %>%
    write_csv("readsout.csv", col_names = F)

# stop(" *** stopped manually *** ") ##########################################
