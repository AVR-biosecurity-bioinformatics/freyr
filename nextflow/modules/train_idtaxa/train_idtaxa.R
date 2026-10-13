#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

ref_fasta <- args$ref_fasta
fasta_type <- args$fasta_type

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "ape",
    "Biostrings",
    "DECIPHER",
    "dplyr",
    "magrittr",
    "purrr",
    "stringr",
    "tibble",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

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

    db <- get_ncbi_taxonomy()
    get_lineage <- TRUE

} else {
    stop("'fasta_type' (params.ncbi_idtaxa) must be 'formatted' or 'unformatted'")
}





idtaxa_model <- train_idtaxa(
    x = hierarchy, 
    max_group_size = 10, 
    max_iterations = 3, 
    allow_group_removal = TRUE, 
    orient = FALSE, 
    get_lineage = get_lineage, 
    db = db, 
    quiet = FALSE
    )

# make filename
out_filename <- stringr::str_remove(basename(ref_fasta), "\\..+$")

# Write out idtaxa objects
saveRDS(idtaxa_model, paste0(out_filename,"_idtaxa_db.rds"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})