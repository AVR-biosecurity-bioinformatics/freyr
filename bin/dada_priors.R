#!/usr/bin/env Rscript

# check variables defined

### run R code
if (direction == "forward") { # recode read direction as "F" or "R"
    direction_short <- "F"
} else if (direction == "reverse") {
    direction_short <- "R"
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

print(priors_list)

priors_tibble <- tibble() # new tibble
for (i in 1:length(priors_list)) { # loop through .rds files, adding distinct sequences to tibble
    seq <- readRDS(priors_list[i])$sequence %>% as_tibble_col(column_name = "sequence") %>% distinct()
    priors_tibble <- rbind(priors_tibble, seq)
}

print(priors_tibble)

# keep only sequences that appear more than once (ie. are in more than one sample)
priors <- priors_tibble %>% 
            group_by(sequence) %>% 
            summarise(n = n()) %>% 
            filter(n>1) %>%
            pull(sequence)

print(priors)

saveRDS(priors, paste0(fcid,"_",pcr_primers,"_priors",direction_short,".rds"))

# stop(" *** stopped manually *** ") ##########################################

