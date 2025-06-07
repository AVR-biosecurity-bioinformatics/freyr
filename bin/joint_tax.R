#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "purrr",
    "readr",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "fcid",
    "pcr_primers",
    "target_gene",
    "idtaxa_list",
    "blast_list"
)
lapply(nf_vars, nf_var_check)

## check and define variables 
target_gene <-          parse_nf_var_repeat(target_gene)

propagate_tax <-        TRUE

### run R code      
idtaxa_output <- 
    idtaxa_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) %>% # read in seqtabs and store as list of tibbles
    purrr::keep(., ~ nrow(.x) > 0) %>% # remove tibbles with 0 rows
    dplyr::bind_rows()

blast_output <- 
    blast_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) %>% # read in seqtabs and store as list of tibbles
    purrr::keep(., ~ nrow(.x) > 0) %>% # remove tibbles with 0 rows
    dplyr::bind_rows()

# Normalise the delimiter in the blast spp output to what is used by idtaxa
if(any(stringr::str_detect(idtaxa_output$Species, "_"), na.rm = TRUE) & !any(stringr::str_detect(idtaxa_output$Species, "_"), na.rm = TRUE)){
  blast_output$blast_spp <- stringr::str_replace_all(blast_output$blast_spp, " ", "_")
} else if(any(str_detect(idtaxa_output$Species, " "), na.rm = TRUE) & !any(stringr::str_detect(idtaxa_output$Species, "_"), na.rm = TRUE)){
  blast_output$blast_spp <- stringr::str_replace_all(blast_output$blast_spp, "_", " ")
}

# check that BLAST ASVs are those in IDTAXA ASVs
if(!all(blast_output$seq_name %in% idtaxa_output$seq_name)){
    stop("ASV names in BLAST output don't match those in IDTAXA output")
}

# try to merge IDTAXA and BLAST assignments
# NOTE: BLAST assignment will only be used if it matches IDTAXA at Genus rank and the IDTAXA Species assignment is NA
if(nrow(idtaxa_output) > 0 & nrow(blast_output) > 0){
    # join idtaxa and blast output together
    joint_assignment <- 
        dplyr::full_join(idtaxa_output, blast_output, by = dplyr::join_by(seq_name)) %>%
        dplyr::mutate(
            # if idtaxa hasn't assigned to Genus level, then set blast assignment to NA
            blast_genus = dplyr::if_else(is.na(Genus), NA, blast_genus),
            blast_spp = dplyr::if_else(is.na(Genus), NA, blast_spp),
            # if blast genus doesn't match idtaxa genus, set former (and blast_spp) to NA
            blast_genus = dplyr::if_else(Genus != blast_genus, NA, blast_genus),
            blast_spp = dplyr::if_else(Genus != blast_genus, NA, blast_spp)
        ) %>%
        # work out the matching between idtaxa and blast assignments at genus and species ranks
        dplyr::mutate(
            species_match = dplyr::case_when(
                is.na(Species) & is.na(blast_spp) ~ "both_NA",
                is.na(Species) & (!is.na(blast_spp)) ~ "idtaxa_NA",
                (!is.na(Species)) & is.na(blast_spp) ~ "blast_NA",
                Species == blast_spp ~ "both_species_match",
                !is.na(Species) & !is.na(blast_spp) ~ "both_species_mismatch",
            )
        ) %>%
        # replace species assignment with BLAST assignment if species_match == "idtaxa_NA"
        dplyr::mutate(Species = dplyr::if_else(species_match == "idtaxa_NA", blast_spp, Species)) %>%
        # keep only taxon rank columns + seq_name
        dplyr::select(seq_name, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
        # replace "NA" taxa with the lowest assigned rank in the "G__Genus" format
        dplyr::mutate(
            # unsure Root is always "Root"
            Root = "Root",
            # determine the lowest assigned rank, choosing the first true case
            lowest_assignment = dplyr::case_when(
                is.na(Kingdom)  ~ paste0("R__",Root),
                is.na(Phylum)   ~ paste0("K__",Kingdom),
                is.na(Class)    ~ paste0("P__",Phylum),
                is.na(Order)    ~ paste0("C__",Class),
                is.na(Family)   ~ paste0("O__",Order),
                is.na(Genus)    ~ paste0("F__",Family),
                is.na(Species)  ~ paste0("G__",Genus),
                .default        = NA
            ),
            # replace NA ranks with the "G__Genus" format ('propagate')
            dplyr::across(
                Kingdom:Species, 
                ~dplyr::case_when(
                    is.na(.) & !is.na(lowest_assignment) ~ lowest_assignment, 
                    .default = .
                )
            )
        ) %>%
        dplyr::select(-lowest_assignment)
} else {
    stop("Either idtaxa_output or blast_output is empty")
}

# another check
if(!all(idtaxa_output$seq_name %in% joint_assignment$seq_name)){
    stop("Not all input ASVs are in the joint assignment output")
}

# Write tibble
readr::write_csv(joint_assignment, paste0(fcid,"_",pcr_primers,"_joint.csv"))

