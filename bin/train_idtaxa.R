#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "ape",
    "Biostrings",
    "DECIPHER",
    "dplyr",
    "magrittr",
    "purrr",
    "stringr",
    "taxreturn",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "pcr_primers",
    "ref_fasta",
    "fasta_type"
)
lapply(nf_vars, nf_var_check)

## check and define variables 

### run R code

# convert fasta to DNABin
hierarchy <- ape::read.FASTA(
    file = ref_fasta,
    type = "DNA"
    )

# formatted vs unformatted fasta input
if ( fasta_type == "formatted" ) {

    db <- NULL
    get_lineage <- FALSE

} else if ( fasta_type == "unformatted" ) {

    db <- taxreturn::get_ncbi_taxonomy()
    get_lineage <- TRUE

} else {
    stop("'fasta_type' (params.ncbi_idtaxa) must be 'formatted' or 'unformatted'")
}





idtaxa_model <- taxreturn::train_idtaxa(
    x = hierarchy, 
    max_group_size = 10, 
    max_iterations = 3, 
    allow_group_removal = TRUE, 
    orient = FALSE, 
    get_lineage = get_lineage, 
    db = db, 
    quiet = FALSE
    )


# Write out idtaxa objects
saveRDS(idtaxa_model, paste0(pcr_primers,"_idtaxa_db.rds"))

# stop(" *** stopped manually *** ") ##########################################
