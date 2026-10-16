#!/usr/bin/env Rscript
.libPaths(.Library)

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

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "phyloseq",
    "readr",
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

# Combine DNAStringSet objects, retaining the first occurrence of each ASV name.
seqs_uf <- unlist(
    lapply(unfiltered_seqs_list, as.character),
    use.names = TRUE
)
seqs_uf <- Biostrings::DNAStringSet(
    seqs_uf[!duplicated(names(seqs_uf))]
)
write_fasta(seqs_uf, "asvs_unfiltered.fasta")


## merge unfiltered phyloseq objects
ps_u <- merge_phyloseq_new(ps_unfiltered)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
expected_ranks <- c(
    "Root",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
)

if (!identical(phyloseq::rank_names(ps_u), expected_ranks)) {
    warning(
        "Taxonomy ranks do not meet staging database expectations.\n",
        "Expected: ", paste(expected_ranks, collapse = ", "), "\n",
        "Observed: ", paste(phyloseq::rank_names(ps_u), collapse = ", ")
    )
}

## output merged data tables 
data.table::fwrite(
    melt_phyloseq(ps_u),
    "raw_unfiltered.csv"
)

# Export species level summary of unfiltered results
summarise_phyloseq(ps_u) %>%
    readr::write_csv(., paste0("summary_unfiltered.csv"))

## output phyloseq and component data

# Write out whole phyloseq object
saveRDS(ps_u, paste0("ps_unfiltered.rds"))

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_primers") %>%
    readr::write_csv("seqtab_unfiltered.csv")

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    ) %>%
    readr::write_csv("taxtab_unfiltered.csv")

# save samplesheet
phyloseq::sample_data(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble() %>%
    readr::write_csv("samdf_unfiltered.csv")

# Clean up
rm(ps_unfiltered)
rm(seqs_uf)
gc()

## read tracking output from unfiltered phyloseq
rank_cols <- colnames(phyloseq::tax_table(ps_u)) # only retrieve this once

# try to save memory
ps_u_summary <- summarise_phyloseq(ps_u)
ps_u_samplenames <- sample_names(ps_u)
ps_u_sampledata <- 
    phyloseq::sample_data(ps_u) %>%
    as("data.frame") %>%
    dplyr::select(sample_primers, read_group, primers)
rm(ps_u)
gc()
message("Pre-readcount files created!")

ps_u_summary %>% 
    dplyr::mutate(
        dplyr::across(
            tidyselect::all_of(ps_u_samplenames), 
            ~ dplyr::na_if(., 0)
        )
    ) %>%
    tidyr::pivot_longer(cols = ps_u_samplenames, names_to="sample_primers", values_to = "Abundance", values_drop_na = TRUE) %>%
    dplyr::left_join(ps_u_sampledata, by = "sample_primers") %>%
    dplyr::mutate(
        dplyr::across(
            tidyselect::any_of(rank_cols), 
            ~ ifelse(!stringr::str_detect(.x, "__"), Abundance, NA_integer_)
        )
    ) %>% 
    dplyr::group_by(sample_primers, read_group, primers) %>%
    dplyr::summarise(
        dplyr::across(
            tidyselect::any_of(rank_cols),
            ~ sum(., na.rm = TRUE),
            .names = "classified_{.col}"
        )
    ) %>% 
    tidyr::pivot_longer(cols = tidyselect::starts_with("classified_"), names_to = "stage", values_to = "pairs") %>% 
    dplyr::mutate(stage = stringr::str_to_lower(stage)) %>%
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>%
    readr::write_csv("ps_u_readsout.csv")

rm(ps_u_summary)
rm(ps_u_samplenames)
rm(ps_u_sampledata)
gc()

### filtered ----------------------------------------------------------------------------------------------------------


ps_filtered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_filtered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_filtered <- lapply(ps_filtered, readRDS) # read in phyloseq objects and store as list

filtered_seqs_list <- 
    filtered_fastas %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)


# Combine DNAStringSet objects, retaining the first occurrence of each ASV name.
seqs_f <- unlist(
    lapply(unfiltered_seqs_list, as.character),
    use.names = TRUE
)
seqs_f <- Biostrings::DNAStringSet(
    seqs_f[!duplicated(names(seqs_f))]
)
write_fasta(seqs_f, "asvs_filtered.fasta")

## merge filtered phyloseq objects
ps_f <- merge_phyloseq_new(ps_filtered)

if (!identical(phyloseq::rank_names(ps_f), expected_ranks)) {
    warning(
        "Taxonomy ranks do not meet staging database expectations.\n",
        "Expected: ", paste(expected_ranks, collapse = ", "), "\n",
        "Observed: ", paste(phyloseq::rank_names(ps_f), collapse = ", ")
    )
}

## output merged raw data tables - NOTE: This is memory intensive
data.table::fwrite(
    melt_phyloseq(ps_f),
    "raw_filtered.csv"
)
# Export species level summary of filtered results
data.table::fwrite(
    summarise_phyloseq(ps_f),
    "summary_filtered.csv"
)

## output phyloseq and component data

# Write out phyloseq
saveRDS(ps_f, paste0("ps_filtered.rds"))

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
seqtab_out_f <- phyloseq::otu_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_primers") %>%
    write_csv("seqtab_filtered.csv")

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out_f <- phyloseq::tax_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    ) %>%
    write_csv("taxtab_filtered.csv")

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out_f) == c("seq_name", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: seq_name, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
phyloseq::sample_data(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble() %>%
    write_csv("samdf_filtered.csv")

## read tracking output from filtered phyloseq
phyloseq::sample_sums(ps_f) %>%
    tibble::enframe(name = "sample_primers", value = "pairs")%>%
    dplyr::left_join(phyloseq::sample_data(ps_f) %>%
                     as("data.frame") %>%
                     dplyr::select(sample_primers, read_group, primers)) %>% 
    dplyr::mutate(stage = "filter_sample_taxon") %>% 
    dplyr::select(stage, sample_primers, read_group, primers, pairs) %>% 
    readr::write_csv("ps_f_readsout.csv")

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})