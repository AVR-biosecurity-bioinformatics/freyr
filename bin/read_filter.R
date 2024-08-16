#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dada2",
    "dplyr",
    "readr",
    "tibble",
    "stringr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "reads_paths",
    "read_min_length",
    "read_max_length",
    "read_max_ee",
    "read_trunc_length",
    "read_trim_left",
    "read_trim_right",
    "sample_id",
    "target_gene",
    "pcr_primers",
    "fcid",
    "seq_type",
    "paired"
)
lapply(nf_vars, nf_var_check)

### process variables 
# split reads_paths into individual file paths (or keep if single reads)
if ( paired == "true" ) {
    reads_paths_vec <- reads_paths %>% stringr::str_split_1(";")
    fwd_reads <- reads_paths_vec[1]
    rev_reads <- reads_paths_vec[2]
} else if ( paired == "false" ) {
   single_reads <- reads_paths
} else {
    stop ( "'paired' must be 'true' or 'false'!" )
}

### run R code

if ( paired == "true" & seq_type == "illumina" ) {

    # filter and trim paired-end reads
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

} else if ( paired == "false" & seq_type == "nanopore" ) {

    # filter and trim single-end reads
    res <- dada2::filterAndTrim(
        fwd = single_reads, 
        filt = paste0(sample_id,"_",target_gene,"_",pcr_primers,"_filter_R0.fastq.gz"), # gets saved to working dir
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

} else {
    stop ( "Currently unsupported sequencing type -- check samplesheet" )
}

## extract output read counts and save to file
reads_out <- res %>%
    tibble::as_tibble() %>%
    dplyr::pull(reads.out)

out_vector <- c("read_filter", sample_id, fcid, pcr_primers, reads_out, reads_out) 

rbind(out_vector) %>%
    tibble::as_tibble() %>%
    readr::write_csv(paste0("read_filter_",sample_id,"_",pcr_primers,"_readsout.csv"), col_names = F)

# stop(" *** stopped manually *** ") ##########################################
