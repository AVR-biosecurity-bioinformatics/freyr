#!/usr/bin/env Rscript

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
    filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    write_csv(., paste0("raw_unfiltered.csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps_u) %>% 
    filter(Abundance > 0) %>%
    left_join(
        refseq(ps_u) %>% as.character() %>% enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, rank_names(ps_u), sample_id, Abundance ) %>%
    pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    write_csv(., paste0("summary_unfiltered.csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps_u)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_unfiltered.fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phy_tree(ps_u, errorIfNULL = FALSE))){
    #Output newick tree
    write.tree(phy_tree(ps_u), file = paste0("tree_unfiltered.nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out_u <- phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out_u <- phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    as_tibble(rownames = "OTU") %>%
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
    as_tibble()

# Write out
write_csv(seqtab_out_u, paste0("seqtab_unfiltered.csv"))
write_csv(taxtab_out_u, paste0("taxtab_unfiltered.csv"))
write_csv(samdf_out_u, paste0("samdf_unfiltered.csv"))
saveRDS(ps_u, paste0("ps_unfiltered.rds"))

## read tracking output from unfiltered phyloseq
phyloseq::psmelt(ps_u) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(sample_id, fcid, Abundance, any_of(colnames(phyloseq::tax_table(ps_u)))) %>%
    pivot_longer(
        cols=any_of(colnames(phyloseq::tax_table(ps_u))), 
        names_to = "rank",
        values_to="name"
        ) %>%
    filter(!is.na(name))%>%
    dplyr::group_by(sample_id, rank) %>%
    summarise(Abundance = sum(Abundance)) %>%
    pivot_wider(
        names_from="rank",
        values_from="Abundance"
        )%>%
    rename_with(~stringr::str_to_lower(.), everything()) %>%
    rename_with(~stringr::str_c("classified_", .), -sample_id) %>% 
    write_csv("ps_u_readsout.csv")


### filtered
## merge filtered phyloseq objects
ps_f <- merge_phyloseq_new(ps_filtered)

## output merged data tables
phyloseq::psmelt(ps_f) %>% 
    filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    write_csv(., paste0("raw_filtered.csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps_f) %>% 
    filter(Abundance > 0) %>%
    left_join(
        refseq(ps_f) %>% as.character() %>% enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, rank_names(ps_f), sample_id, Abundance ) %>%
    pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    write_csv(., paste0("summary_filtered.csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps_f)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_filtered.fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phy_tree(ps_f, errorIfNULL = FALSE))){
    #Output newick tree
    write.tree(phy_tree(ps_f), file = paste0("tree_filtered.nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out_f <- phyloseq::otu_table(ps_f) %>%
    as("matrix") %>%
    as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out_f <- phyloseq::tax_table(ps_f) %>%
    as("matrix") %>%
    as_tibble(rownames = "OTU") %>%
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
    as_tibble()

# Write out
write_csv(seqtab_out_f, paste0("seqtab_filtered.csv"))
write_csv(taxtab_out_f, paste0("taxtab_filtered.csv"))
write_csv(samdf_out_f, paste0("samdf_filtered.csv"))
saveRDS(ps_f, paste0("ps_filtered.rds"))




# stop(" *** stopped manually *** ") ##########################################
