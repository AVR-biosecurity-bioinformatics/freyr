#!/usr/bin/env Rscript


## check and define variables

#### INCOMPLETE
fcid_list <- # convert Groovy to R list format
    stringr::str_extract_all(fcid, pattern = "[^\\s,\\[\\]]+") %>% unlist()

pcr_primers_list <- # convert Groovy to R list format
    stringr::str_extract_all(pcr_primers, pattern = "[^\\s,\\[\\]]+") %>% unlist()

meta_list <-  # convert Groovy to R list format
    stringr::str_extract_all(meta, pattern = "[^\\[\\]]+") %>% 
    unlist() %>% str_subset(pattern = "^[^,]") # remove elements that start with ","
meta_list2 <- meta_list %>% stringr::str_split(pattern = ", ") # nested list; elements are samples x pcr_primers

taxtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(taxtab, pattern = "[^\\s,\\[\\]]+") %>% unlist()

taxtab_tibble <- tibble() # read in taxtabs
for (i in 1:length(taxtab_list)) { # loop through .rds files, adding distinct sequences to tibble
    taxtab_i <- readRDS(taxtab_list[i])
    taxtab_tibble <- rbind(taxtab_tibble, taxtab_i)

}
#### INCOMPLETE



### run R code