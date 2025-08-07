#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

ps_unfiltered               <- args$ps_unfiltered
ps_filtered                 <- args$ps_filtered
unfiltered_fastas           <- args$unfiltered_fastas
filtered_fastas             <- args$filtered_fastas

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "phyloseq",
    "readr",
    "seqateurs",
    "stringr",
    "tibble",
    "data.table",
    "tidyr",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

## check and define variables
ps_unfiltered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_unfiltered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_unfiltered <- lapply(ps_unfiltered, readRDS) # read in phyloseq objects and store as list

# import fasta files
unfiltered_seqs_list <- 
    unfiltered_fastas %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

### run R code

### unfiltered ---------------------------------------------------------------------------------------------

# combine fasta DSS objects, removing redundant sequences
seqs_uf <- 
    unfiltered_seqs_list %>%
    lapply(., as.character) %>%
    unlist(use.names = T) %>% 
    tibble::enframe() %>%
    dplyr::group_by(name, value) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

# save combined fasta of filtered sequences
write_fasta(seqs_uf, "asvs_unfiltered.fasta")

gc()

## merge unfiltered phyloseq objects
ps_u <- merge_phyloseq_new(ps_unfiltered)

gc()

## output merged data tables 
melt_phyloseq(ps_u) %>% 
  readr::write_csv(., paste0("raw_unfiltered.csv"))

gc()

# Export species level summary of unfiltered results
summarise_phyloseq(ps_u) %>%
    readr::write_csv(., paste0("summary_unfiltered.csv"))

gc()

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
seqtab_out_u <- phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_primers")

gc()

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out_u <- phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    )

gc()

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out_u) == c("seq_name", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: seq_name, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
samdf_out_u <- phyloseq::sample_data(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble()

# Write out
readr::write_csv(seqtab_out_u, paste0("seqtab_unfiltered.csv"))
readr::write_csv(taxtab_out_u, paste0("taxtab_unfiltered.csv"))
readr::write_csv(samdf_out_u, paste0("samdf_unfiltered.csv"))
saveRDS(ps_u, paste0("ps_unfiltered.rds"))

gc()

## read tracking output from unfiltered phyloseq
rank_cols <- colnames(phyloseq::tax_table(ps_u)) # only retrieve this once

summarise_phyloseq(ps_u) %>%
    tidyr::pivot_longer(cols = sample_names(ps_u), names_to="sample_primers", values_to = "Abundance")%>%
    dplyr::left_join(
        phyloseq::sample_data(ps_u) %>%
            as("data.frame") %>%
            dplyr::select(sample_primers, read_group, primers)
        ) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(rank_cols), 
                                ~ ifelse(!stringr::str_detect(.x, "__"), Abundance, NA_integer_))) %>% 
    dplyr::group_by(sample_primers, read_group, primers) %>%
    dplyr::summarise(dplyr::across(tidyselect::any_of(rank_cols),
                                    ~ sum(., na.rm = TRUE), .names = "classified_{.col}")) %>% 
    tidyr::pivot_longer(cols = tidyselect::starts_with("classified_"), names_to = "stage", values_to = "pairs") %>% 
    mutate(stage = stringr::str_to_lower(stage)) %>%
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>%
    readr::write_csv("ps_u_readsout.csv")

rm(ps_u)
rm(ps_unfiltered)
rm(seqs_uf)
rm(seqtab_out_u)
rm(taxtab_out_u)
rm(samdf_out_u)
gc()

### filtered ----------------------------------------------------------------------------------------------------------


ps_filtered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_filtered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_filtered <- lapply(ps_filtered, readRDS) # read in phyloseq objects and store as list

filtered_seqs_list <- 
    filtered_fastas %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

# combine fasta DSS objects, removing redundant sequences
seqs_f <- 
    filtered_seqs_list %>%
    lapply(., as.character) %>%
    unlist(use.names = T) %>% 
    tibble::enframe() %>%
    dplyr::group_by(name, value) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

# save combined fasta of filtered sequences
write_fasta(seqs_f, "asvs_filtered.fasta")

gc()

## merge filtered phyloseq objects
ps_f <- merge_phyloseq_new(ps_filtered)

gc()

## output merged raw data tables - NOTE: This is memory intensive
melt_phyloseq(ps_f) %>% 
    readr::write_csv(., paste0("raw_filtered.csv"))

gc()

# Export species level summary of filtered results
summarise_phyloseq(ps_f) %>%
    readr::write_csv(., paste0("summary_filtered.csv"))

gc()

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
seqtab_out_f <- phyloseq::otu_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_primers")

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out_f <- phyloseq::tax_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    )

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out_f) == c("seq_name", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: seq_name, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
samdf_out_f <- phyloseq::sample_data(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble()

# Write out
readr::write_csv(seqtab_out_f, paste0("seqtab_filtered.csv"))
readr::write_csv(taxtab_out_f, paste0("taxtab_filtered.csv"))
readr::write_csv(samdf_out_f, paste0("samdf_filtered.csv"))
saveRDS(ps_f, paste0("ps_filtered.rds"))

gc()

## read tracking output from filtered phyloseq
phyloseq::sample_sums(ps_f) %>%
    tibble::enframe(name = "sample_primers", value = "pairs")%>%
    dplyr::left_join(phyloseq::sample_data(ps_f) %>%
                     as("data.frame") %>%
                     dplyr::select(sample_primers, read_group, primers)) %>% 
    dplyr::mutate(stage = "filter_sample_taxon") %>% 
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>% 
    readr::write_csv("ps_f_readsout.csv")
    

gc()

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})