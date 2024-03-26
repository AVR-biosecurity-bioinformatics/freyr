#!/usr/bin/env Rscript

# check variables defined


### run R code

# # using piperline function
# step_filter_reads(
#     sample_id = ,
#     input_dir = ,
#     output_dir = ,
#     min_length = ,
#     max_length = ,


# )

# using dada2::filterAndTrim directly
dada2::filterAndTrim(
    fwd = fwd_reads, 
    filt = "R1_out.fastq.gz",
    rev = rev_reads, 
    filt.rev = "R2_out.fastq.gz",
    minLen = read_min_length, 
    maxLen = read_max_length, 
    maxEE = read_max_ee, 
    truncLen = read_trunc_length,
    trimLeft = read_trim_left, 
    trimRight = read_trim_right, 
    rm.phix = TRUE, 
    multithread = FALSE, 
    compress = TRUE, 
    verbose = !quiet
)