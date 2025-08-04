#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads_paths                 <- args$reads_paths
read_min_length             <- args$read_min_length
read_max_length             <- args$read_max_length
read_max_ee                 <- args$read_max_ee
read_trunc_length           <- args$read_trunc_length
read_trim_left              <- args$read_trim_left
read_trim_right             <- args$read_trim_right
sample_primers              <- args$sample_primers
locus                       <- args$locus
primers                     <- args$primers
read_group                  <- args$read_group
seq_type                    <- args$seq_type
paired                      <- args$paired

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dada2",
    "dplyr",
    "readr",
    "tibble",
    "stringr",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

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

    ## check if input files have reads in them
    if (file.info(fwd_reads)$size > 0 && file.info(rev_reads)$size > 0){
        # filter and trim paired-end reads
        set.seed(1); res <- dada2::filterAndTrim(
            fwd = fwd_reads, 
            filt = paste0(sample_primers,"_",locus,"_",primers,"_filter_R1.fastq.gz"), # gets saved to working dir
            rev = rev_reads, 
            filt.rev = paste0(sample_primers,"_",locus,"_",primers,"_filter_R2.fastq.gz"), # gets saved to working dir
            minLen = as.numeric(read_min_length), 
            maxLen = as.numeric(read_max_length), 
            maxEE = as.numeric(read_max_ee), 
            truncLen = as.numeric(read_trunc_length),
            trimLeft = as.numeric(read_trim_left), 
            trimRight = as.numeric(read_trim_right), 
            rm.phix = TRUE, 
            multithread = FALSE, 
            compress = TRUE, 
            verbose = FALSE,
            n = 1e5
        )
    } else {
        message("Input file size is 0, skipping read trimming")
    }

    # create output files if they don't exist
    if (!file.exists(paste0(sample_primers,"_",locus,"_",primers,"_filter_R1.fastq.gz"))){
        file.create(paste0(sample_primers,"_",locus,"_",primers,"_filter_R1.fastq.gz"))
        empty_forward <- TRUE
    } else {
        empty_forward <- FALSE
    }
    if (!file.exists(paste0(sample_primers,"_",locus,"_",primers,"_filter_R2.fastq.gz"))){
        file.create(paste0(sample_primers,"_",locus,"_",primers,"_filter_R2.fastq.gz"))
        empty_reverse <- TRUE
    } else {
        empty_reverse <- FALSE
    }
    if (empty_forward & empty_reverse ){
        empty_reads <- TRUE
    } else if ( !empty_forward & !empty_reverse ) {
        empty_reads <- FALSE
    } else {
        stop(paste0("Forward reads empty is ",empty_forward," while reverse read empty is ",empty_reverse))
    }

} else if ( paired == "false" & seq_type == "nanopore" ) {
    
    ## check if input file has reads in it
    if (file.info(single_reads)$size > 0){
        # filter and trim single-end reads
        set.seed(1); res <- dada2::filterAndTrim(
            fwd = single_reads, 
            filt = paste0(sample_primers,"_",locus,"_",primers,"_filter_R0.fastq.gz"), # gets saved to working dir
            minLen = as.numeric(read_min_length), 
            maxLen = as.numeric(read_max_length), 
            maxEE = as.numeric(read_max_ee), 
            truncLen = as.numeric(read_trunc_length),
            trimLeft = as.numeric(read_trim_left), 
            trimRight = as.numeric(read_trim_right), 
            rm.phix = TRUE, 
            multithread = FALSE, 
            compress = TRUE, 
            verbose = FALSE,
            n = 1e4
        )
    } else {
        message("Input file size is 0, skipping read trimming")
    }

    # create output file if they don't exist
    if (!file.exists(paste0(sample_primers,"_",locus,"_",primers,"_filter_R0.fastq.gz"))){
        file.create(paste0(sample_primers,"_",locus,"_",primers,"_filter_R0.fastq.gz"))
        empty_reads <- TRUE
    } else {
        empty_reads <- FALSE
    }

} else {
    stop ( "Currently unsupported sequencing type -- check samplesheet" )
}


## extract output read counts and save to file
if ( !empty_reads ){
    reads_out <- res %>%
        tibble::as_tibble() %>%
        dplyr::pull(reads.out)
} else {
    reads_out <- 0
}
out_vector <- c("read_filter", sample_primers, read_group, primers, reads_out, reads_out) 

rbind(out_vector) %>%
    tibble::as_tibble() %>%
    readr::write_csv(paste0("read_filter_",sample_primers,"_",primers,"_readsout.csv"), col_names = F)


# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})