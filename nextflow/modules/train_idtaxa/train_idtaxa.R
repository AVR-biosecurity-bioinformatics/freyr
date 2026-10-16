#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

ref_fasta          <- args$ref_fasta
max_group_size     <- as.integer(args$max_group_size)
max_iterations     <- as.integer(args$max_iterations)
allow_group_removal <- as.logical(args$allow_group_removal)
drop_malformed     <- as.logical(args$drop_malformed)

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

# Train IDTAXA model
idtaxa_model <- train_idtaxa(
    x = hierarchy, 
    max_group_size = max_group_size, 
    max_iterations = max_iterations, 
    allow_group_removal = allow_group_removal, 
    drop_malformed = drop_malformed,
    orient = FALSE, 
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