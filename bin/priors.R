#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
direction                   <- args$direction
read_group                  <- args$read_group
priors                      <- args$priors
prior_min_samples           <- args$prior_min_samples

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### run R code
if (direction == "forward") { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if (direction == "reverse") {
    direction_short <- "R"
} else if ( direction == "single" ) {
    direction_short <- "S"
} else {
    stop(" Input reads direction needs to be 'forward' or 'reverse'! ")
}

prior_min_samples <- as.numeric(prior_min_samples)

## parse priors list
priors_list <- # convert input priors .rds file list from Groovy format to R format
    stringr::str_extract_all(
        priors, 
        pattern = "\\S+?\\.rds" 
        ) %>% 
    unlist() 

priors_tibble <- tibble::tibble() # new tibble
for (i in 1:length(priors_list)) { # loop through .rds files, adding distinct sequences to tibble
    seq <- readRDS(priors_list[i])$sequence %>% tibble::as_tibble_col(column_name = "sequence") %>% dplyr::distinct()
    priors_tibble <- rbind(priors_tibble, seq)
}
if (nrow(priors_tibble) > 0){
    # keep only sequences that appear in at least --priors_min_samples samples
    priors <- priors_tibble %>% 
        dplyr::group_by(sequence) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::filter(n >= prior_min_samples) %>%
        dplyr::pull(sequence)

    # check there are sequences remaining after filtering
    if ( identical(priors, character(0)) ){
        stop(paste0("\n***\nThere were ",nrow(priors_tibble)," input prior sequences for read group '",read_group,"' and primers '", primers, "', but all of them were removed during filtering.\nConsider decreasing the value of --prior_min_samples, or disabling --high_sensitivity, or checking that enough reads are passing through the read filtering steps.\n***"))
    }

} else {
    priors <- NULL
}

saveRDS(priors, paste0(read_group,"_",primers,"_priors",direction_short,".rds"))

# stop(" *** stopped manually *** ") ##########################################

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})