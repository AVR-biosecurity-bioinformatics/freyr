#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "purrr",
    "readr",
    "stringr",
    "tidyr",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "samplesheet",
    "loci_params"
)
lapply(nf_vars, nf_var_check)

### run code

samplesheet_df <- readr::read_csv(paste0(projectDir,"/",samplesheet), show_col_types = F) # read in samplesheet

loci_params_df <- readr::read_csv(paste0(projectDir,"/",loci_params), show_col_types = F) # read in loci parameters


## validation of samplesheet and loci_params content

# check that sample_name is a subset of sample_id

# check that sample_id contains only unique values
if (any(duplicated(samplesheet_df$sample_id))) { stop("INPUT ERROR: All sample IDs in samplesheet must be unique!") } 

# check that all primer names and sequences are identical for each sample
if (!any(duplicated(samplesheet_df$pcr_primers))) { stop("INPUT ERROR: 'pcr_primer' varies between samples!") } 
if (!any(duplicated(samplesheet_df$for_primer_seq))) { stop("INPUT ERROR: 'for_primer_seq' varies between samples!") } 
if (!any(duplicated(samplesheet_df$rev_primer_seq))) { stop("INPUT ERROR: 'rev_primer_seq' varies between samples!") } 

# check that same number of target genes and primer names are defined as primer sequences
if (
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) ||
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("INPUT ERROR: Mismatch between number of supplied target genes and primer sequences!") }
if (
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) ||
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("INPUT ERROR: Mismatch between number of supplied primer pair names and primer sequences!") }

### split samplesheets into loci


stop(" *** stopped manually *** ") ##########################################
