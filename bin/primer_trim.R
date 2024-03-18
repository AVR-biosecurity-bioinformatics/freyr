#!/usr/bin/env Rscript

### this runs step_primer_trim() on samples per primer pair
## inputs needs to be a single read pair

### inputs:
## sample_id
## input_dir
## output_dir
## qc_dir
## for_primer_seq
## rev_primer_seq
## pcr_primers

step_primer_trim(
    sample_id = ..1, 
    for_primer_seq=..2, 
    rev_primer_seq=..3, 
    pcr_primers = ..4,
    input_dir = paste0("data/",..5), 
    output_dir =  paste0("data/",..5,"/01_trimmed"),
    qc_dir=paste0("output/logs/",..5),
    max_mismatch=..6,
    quiet = FALSE)