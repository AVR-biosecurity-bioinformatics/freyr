#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "direction",
    "read_group",
    "primers",
    "priors"
)
lapply(nf_vars, nf_var_check)

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

# keep only sequences that appear more than once (ie. are in more than one sample)
priors <- priors_tibble %>% 
            dplyr::group_by(sequence) %>% 
            dplyr::summarise(n = n()) %>% 
            dplyr::filter(n>1) %>%
            dplyr::pull(sequence)

saveRDS(priors, paste0(read_group,"_",primers,"_priors",direction_short,".rds"))

# stop(" *** stopped manually *** ") ##########################################

