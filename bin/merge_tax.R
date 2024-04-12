#!/usr/bin/env Rscript


## check and define variables
fcid_list <- # convert Groovy to R list format
    stringr::str_extract_all(
            fcid, 
            pattern = "[^\\s,\\[\\]]+"
            ) %>% 
        unlist()

pcr_primers_list <- # convert Groovy to R list format
    stringr::str_extract_all(
            pcr_primers, 
            pattern = "[^\\s,\\[\\]]+"
            ) %>% 
        unlist()

meta_list <-  # convert Groovy to R list format
stringr::str_extract_all(
            meta, 
            pattern = "[^\\[\\]]+" # pull out lists (ignore nesting structure)
            ) %>% 
    unlist() %>% 
    str_subset(pattern = "^[^,]") # remove elements that start with ","


taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(
            taxtab, 
            pattern = "[^\\s,\\[\\]]+"
            ) %>% 
        unlist()

taxtab_tibble <- tibble() # read in taxtabs
for (i in 1:length(taxtab_list)) { # loop through .rds files, adding distinct sequences to tibble
    taxtab_i <- readRDS(taxtab_list[i])
    taxtab_tibble <- rbind(taxtab_tibble, taxtab_i)





### run R code