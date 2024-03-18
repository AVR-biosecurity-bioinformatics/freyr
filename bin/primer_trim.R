#!/usr/bin/env Rscript

### this runs step_primer_trim() on samples per primer pair
## inputs need to be a single read pair, which already have paths from channel
## also metadata needs to be input

### inputs:
## sample_id
## input_dir
## output_dir
## qc_dir
## for_primer_seq
## rev_primer_seq
## pcr_primers

input_list <- c(
    sample_id,
    for_primer_seq,
    rev_primer_seq,
    pcr_primers,
    fcid,
    output_dir,
    qc_dir,
    max_mismatch
)

for (i in input_list) {
    print(i)
}

exit(status = 1)

data_loc_abs <- paste0(projectDir,"/",data_loc)


step_primer_trim(
    sample_id = sample_id, 
    for_primer_seq = for_primer_seq, 
    rev_primer_seq = rev_primer_seq, 
    pcr_primers = pcr_primers,
    input_dir = paste0(projectDir,"/",data_loc,"/",fcid), 
    output_dir =  paste0(projectDir,"/",data_loc,"/",fcid,"/01_trimmed"),
    qc_dir=paste0(projectDir,"/output/logs/",fcid),
    max_mismatch=..6,
    quiet = FALSE
    )

### OR use trim_primers() directly

# check inputs
# check input_dir
# create output_dir if needed
# check qc_dir exists and create if not
# check both primers defined
# define output file names and paths


# trim primers
trim_primers(
    fwd = , # file path for fwd read
    rev = , # file path for rev read
    fwd_out = , # file path and name for fwd outfile
    rev_out= , # file path and name for rev outfile
    for_primer_seq = , # sequence of fwd primer; code has "stringr::str_replace_all(for_primer_seq, "I", "N"),"
    rev_primer_seq = , # sequence of rev primer; dode has "stringr::str_replace_all(rev_primer_seq, "I", "N"),"
    n = 1e6, # number of reads to stream at a time?
    qualityType = "Auto", 
    check_paired = TRUE, 
    id.field = NULL, 
    max_mismatch=0, 
    id.sep = "\\s", 
    compress = TRUE, 
    quiet = FALSE
)