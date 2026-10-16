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
max_group_size <- args$max_group_size
max_iterations <- args$max_iterations
allow_group_removal <- as.logical(args$allow_group_removal)
drop_malformed <- as.logical(args$drop_malformed)

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

# Validate format is correct
if (length(hierarchy) == 0) {
    stop("Reference FASTA contains no sequences.")
}

headers <- names(hierarchy)

if (is.null(headers) || any(headers == "")) {
    stop("All reference sequences must have FASTA headers.")
}

if (any(!grepl(";", headers, fixed = TRUE))) {
    stop("All FASTA headers must contain a semicolon-delimited taxonomy.")
}

if (any(grepl(";;", headers, fixed = TRUE))) {
    stop("Empty taxonomy ranks detected (';;').")
}

na_taxonomy <- grepl("(^|;)NA(;|$)", headers)

if (any(na_taxonomy)) {
    bad_headers <- headers[na_taxonomy]
    max_print <- 20L

    message(
        length(bad_headers),
        " FASTA record(s) contain an NA taxonomy rank:\n",
        paste0("  ", head(bad_headers, max_print), collapse = "\n"),
        if (length(bad_headers) > max_print) {
            paste0(
                "\n  ... and ",
                length(bad_headers) - max_print,
                " additional record(s)"
            )
        } else {
            ""
        }
    )

    if (drop_malformed) {
        hierarchy <- hierarchy[!na_taxonomy]
        headers <- names(hierarchy)

        message(
            "Dropped ",
            sum(na_taxonomy),
            " malformed FASTA record(s); ",
            length(hierarchy),
            " record(s) remain."
        )

        if (length(hierarchy) == 0L) {
            stop(
                "No reference sequences remain after dropping malformed records."
            )
        }
    } else {
        stop(
            "One or more taxonomy ranks are NA. ",
            "Set '--drop_malformed true' to remove these records."
        )
    }
}

# Train IDTAXA model
idtaxa_model <- train_idtaxa(
    seqs = hierarchy, 
    max_group_size = max_iterations, 
    max_iterations = max_iterations, 
    allow_group_removal = allow_group_removal, 
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