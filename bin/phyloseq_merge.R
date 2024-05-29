#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "phyloseq",
    "readr",
    "seqateurs",
    "stringr",
    "tibble",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "ps_unfiltered",
    "ps_filtered"
)
lapply(nf_vars, nf_var_check)

## check and define variables
ps_unfiltered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_unfiltered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_unfiltered <- lapply(ps_unfiltered, readRDS) # read in phyloseq objects and store as list


ps_filtered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_filtered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_filtered <- lapply(ps_filtered, readRDS) # read in phyloseq objects and store as list

### run R code

### unfiltered 
## merge unfiltered phyloseq objects
ps_u <- merge_phyloseq_new(ps_unfiltered)

## output merged data tables
phyloseq::psmelt(ps_u) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    readr::write_csv(., paste0("raw_unfiltered.csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps_u) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::left_join(
        phyloseq::refseq(ps_u) %>% as.character() %>% tibble::enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, phyloseq::rank_names(ps_u), sample_id, Abundance ) %>%
    tidyr::pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    readr::write_csv(., paste0("summary_unfiltered.csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps_u)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_unfiltered.fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phyloseq::phy_tree(ps_u, errorIfNULL = FALSE))){
    #Output newick tree
    ape::write.tree(phyloseq::phy_tree(ps_u), file = paste0("tree_unfiltered.nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out_u <- phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out_u <- phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    seqateurs::unclassified_to_na(rownames = FALSE)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out_u) == c("OTU", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: OTU, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
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

## read tracking output from unfiltered phyloseq
rank_cols <- colnames(phyloseq::tax_table(ps_u)) # only retrieve this once
rank_cols_lower <- stringr::str_to_lower(rank_cols)

phyloseq::psmelt(ps_u) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(sample_id, fcid, pcr_primers, Abundance, tidyselect::any_of(rank_cols)) %>%
    tidyr::pivot_longer(
        cols=tidyselect::any_of(colnames(phyloseq::tax_table(ps_u))), 
        names_to = "rank",
        values_to="name"
        ) %>%
    dplyr::filter(!is.na(name))%>%
    dplyr::group_by(sample_id, fcid, pcr_primers, rank) %>%
    dplyr::summarise(Abundance = sum(Abundance)) %>%
    tidyr::pivot_wider(
        names_from="rank",
        values_from="Abundance"
        )%>%
    dplyr::rename_with(~stringr::str_to_lower(.), tidyselect::any_of(rank_cols)) %>%
    dplyr::rename_with(~stringr::str_c("classified_", .), tidyselect::any_of(rank_cols_lower)) %>% 
    tidyr::pivot_longer(cols = tidyselect::starts_with("classified_"), names_to = "stage", values_to = "pairs") %>% 
    dplyr::select(stage, sample_id, fcid, pcr_primers, pairs) %>% 
    readr::write_csv("ps_u_readsout.csv")


### filtered
## merge filtered phyloseq objects
ps_f <- merge_phyloseq_new(ps_filtered)

## output merged data tables
phyloseq::psmelt(ps_f) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    readr::write_csv(., paste0("raw_filtered.csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps_f) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::left_join(
        phyloseq::refseq(ps_f) %>% as.character() %>% tibble::enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, phyloseq::rank_names(ps_f), sample_id, Abundance ) %>%
    tidyr::pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    readr::write_csv(., paste0("summary_filtered.csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps_f)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_filtered.fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phyloseq::phy_tree(ps_f, errorIfNULL = FALSE))){
    #Output newick tree
    ape::write.tree(phyloseq::phy_tree(ps_f), file = paste0("tree_filtered.nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out_f <- phyloseq::otu_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out_f <- phyloseq::tax_table(ps_f) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    seqateurs::unclassified_to_na(rownames = FALSE)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out_f) == c("OTU", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: OTU, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
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

## read tracking output from filtered phyloseq
phyloseq::psmelt(ps_f) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(sample_id, fcid, pcr_primers, Abundance) %>%
    dplyr::group_by(sample_id, fcid, pcr_primers) %>%
    dplyr::summarise(pairs = sum(Abundance)) %>% 
    dplyr::mutate(stage = "filter_sample_taxon") %>% 
    dplyr::select(stage, sample_id, fcid, pcr_primers, pairs) %>% 
    readr::write_csv("ps_f_readsout.csv")


# stop(" *** stopped manually *** ") ##########################################
