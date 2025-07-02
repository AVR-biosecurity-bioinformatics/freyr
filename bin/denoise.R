#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads                       <- args$reads
primers                     <- args$primers
sample_primers              <- args$sample_primers
read_group                  <- args$read_group
direction                   <- args$direction
threads                     <- args$cpus
errormodel_file             <- args$errormodel
n_pass                      <- args$n_pass
priors_file                 <- args$priors
dada_band_size              <- args$dada_band_size
dada_homopolymer            <- args$dada_homopolymer

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dada2",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### parse variables
if ( dada_band_size == "null" ){
    dada_band_size <- 16 # default value
} else {
    dada_band_size <- as.numeric(dada_band_size)
}

if ( dada_homopolymer == "null" ){
    dada_homopolymer <- NULL
} else {
    dada_homopolymer <- as.numeric(dada_homopolymer) 
}

errormodel <- readRDS(errormodel_file) # import error model

if ( direction == "forward" ) { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if ( direction == "reverse" ) {
    direction_short <- "R"
} else if ( direction == "single" ) {
    direction_short <- "S"
} else {
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

### run R code

if ( n_pass == "first" && priors_file == "NO_FILE" ) { # first pass condition; no priors
    priors <- character(0)

    # run dada2
    set.seed(1); dada_output <- dada2::dada(
        derep = reads, 
        err = errormodel, 
        multithread = FALSE, 
        priors = priors, 
        selfConsist = FALSE, 
        pool = FALSE, 
        verbose = TRUE,
        BAND_SIZE = dada_band_size,
        HOMOPOLYMER_GAP_PENALTY = dada_homopolymer
    )

    saveRDS(dada_output, paste0(sample_primers,"_dada1",direction_short,".rds"))

} else if ( n_pass == "second" && !priors_file == "NO_FILE" ) { # second pass condition; priors included
    # import priors
    priors <- readRDS(priors_file)

    # run dada2
    set.seed(1); dada_output <- dada2::dada(
        derep = reads, 
        err = errormodel, 
        multithread = FALSE, 
        priors = priors, 
        selfConsist = FALSE, 
        pool = FALSE, 
        verbose = TRUE,
        BAND_SIZE = dada_band_size,
        HOMOPOLYMER_GAP_PENALTY = dada_homopolymer
    )

    saveRDS(dada_output, paste0(sample_primers,"_dada2",direction_short,".rds"))

} else {
    stop(" 'n_pass' variable must be 'first' or 'second', and priors must be 'NO_FILE' or defined! ")
}

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})