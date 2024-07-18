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
    "pcr_primers",
    "ps",
    "target_kingdom",
    "target_phylum",
    "target_class",
    "target_order",
    "target_family",
    "target_genus",
    "target_species",
    "min_sample_reads",
    "min_taxa_reads",
    "min_taxa_ra"
)
lapply(nf_vars, nf_var_check)

## check and define variables
ps <- readRDS(ps)

# convert "NA" strings to true NA
if(target_kingdom == "NA"){ target_kingdom <- NA }
if(target_phylum == "NA"){ target_phylum <- NA }
if(target_class == "NA"){ target_class <- NA }
if(target_order == "NA"){ target_order <- NA }
if(target_family == "NA"){ target_family <- NA }
if(target_genus == "NA"){ target_genus <- NA }
if(target_species == "NA"){ target_species <- NA }

# convert "NA" strings to true NA, or convert number strings to numeric
if(min_sample_reads != "NA"){ min_sample_reads <- as.numeric(min_sample_reads) } else { min_sample_reads <- NA }
if(min_taxa_reads != "NA"){ min_taxa_reads <- as.numeric(min_taxa_reads) } else { min_taxa_reads <- NA }
if(min_taxa_ra != "NA"){ min_taxa_ra <- as.numeric(min_taxa_ra) } else { min_taxa_ra <- NA }

quiet <- FALSE

### run R code

## from step_filter_phyloseq()
# Taxonomic filtering
taxtab <- phyloseq::tax_table(ps) %>%
    as("matrix") %>% 
    as.data.frame()

# Check if any taxonomic filters are enabled
ps0 <- ps

if(any(!sapply(c(target_kingdom, target_phylum, target_class, target_order, target_family, target_genus, target_species), is.na))){

    # Filter kingdom
    if(!is.na(target_kingdom)){
        if (any(stringr::str_detect(taxtab$Kingdom, target_kingdom))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Kingdom",
                value = target_kingdom
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Kingdom", target_kingdom," - Check your target_kingdom parameter"))
        }
    }
    # Filter phylum
    if(!is.na(target_phylum)){
        if (any(stringr::str_detect(taxtab$Phylum, target_phylum))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Phylum",
                value = target_phylum
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Phylum", target_phylum," - Check your target_phylum parameter"))
        }
    }
    # Filter class
    if(!is.na(target_class)){
        if (any(stringr::str_detect(taxtab$Class, target_class))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Class",
                value = target_class
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Class", target_class," - Check your target_class parameter"))
        }
    }
    # Filter order
    if(!is.na(target_order)){
        if (any(stringr::str_detect(taxtab$Order, target_order))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Order",
                value = target_order
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Order", target_order," - Check your target_order parameter"))
        }
    }
    # Filter family
    if(!is.na(target_family)){
        if (any(stringr::str_detect(taxtab$Family, target_family))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Family",
                value = target_family
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Family", target_family," - Check your target_family parameter"))
        }
    }
    # Filter genus
    if(!is.na(target_genus)){
        if (any(stringr::str_detect(taxtab$Genus, target_genus))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Genus",
                value = target_genus
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Genus", target_genus," - Check your target_genus parameter"))
        }
    }
    # Filter Species
    if(!is.na(target_species)){
        if (any(stringr::str_detect(taxtab$Species, target_species))){
        ps0 <- ps0 %>%
            subset_taxa_new(
                rank = "Species",
                value = target_species
            ) %>%
            phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE)
        } else{
        warning(paste0("No ASVs were assigned to the Species", target_species," - Check your target_species parameter"))
        }
    }
} else {
    if (!quiet){message(paste0("No taxonomic filters set - skipping this filter"))}
}

# Remove any taxa under read count or relative abundance thresholds
### TODO: Add comments explaining what is happening here
if(!is.na(min_taxa_reads) & is.na(min_taxa_ra)){
    ps1 <- phyloseq::transform_sample_counts(ps0, function(OTU, ab = min_taxa_reads){ ifelse(OTU <= ab,  0, OTU) })
} else if(is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps1 <- phyloseq::transform_sample_counts(ps0, function(OTU, ab = min_taxa_ra ){ ifelse((OTU / sum(OTU)) <= ab,  0, OTU) })
} else if (!is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps1 <- ps0 %>%
        phyloseq::transform_sample_counts(function(OTU, ab = min_taxa_reads){ ifelse(OTU <= ab,  0, OTU) }) %>%
        phyloseq::transform_sample_counts(function(OTU, ab = min_taxa_ra ){ ifelse((OTU / sum(OTU)) <= ab,  0, OTU) })
} else {
    if (!quiet){message(paste0("No minimum abundance filters set - skipping this filter"))}
    ps1 <- ps0
}

## check if any 


#Remove all samples under the minimum read threshold 
if(min_sample_reads > 0){
    if (all(phyloseq::sample_sums(ps1)<min_sample_reads)) {
        stop(paste0("ERROR: No samples contained reads above the minimum threshold of ", min_sample_reads, " -- consider lowering this value"))
        }
    
    ps2 <- ps1 %>%
        phyloseq::prune_samples(phyloseq::sample_sums(.)>=min_sample_reads, .) %>% 
        phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
} else {
    if (!quiet){message(paste0("No minimum sample reads filter set - skipping this filter"))}
    ps2 <- ps1 %>%
        phyloseq::prune_samples(phyloseq::sample_sums(.)>=0,.) %>%
        phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
}

# Message how many were removed
if(!quiet){message(phyloseq::nsamples(ps) - phyloseq::nsamples(ps2), " Samples and ", phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps2), " ASVs dropped")}

### output filtered results per locus; from step_output_summary()
# Export raw csv
phyloseq::psmelt(ps2) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    readr::write_csv(., paste0("raw_filtered_",pcr_primers,".csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps2) %>% 
    dplyr::filter(Abundance > 0) %>%
    dplyr::left_join(
        phyloseq::refseq(ps2) %>% as.character() %>% tibble::enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, rank_names(ps2), sample_id, Abundance ) %>%
    tidyr::pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    readr::write_csv(., paste0("summary_filtered_",pcr_primers,".csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps2)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_filtered_",pcr_primers,".fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phyloseq::phy_tree(ps2, errorIfNULL = FALSE))){
    #Output newick tree
    ape::write.tree(phyloseq::phy_tree(ps2), file = paste0("tree_filtered",pcr_primers,".nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out <- phyloseq::otu_table(ps2) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out <- phyloseq::tax_table(ps2) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    seqateurs::unclassified_to_na(rownames = FALSE)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out) == c("OTU", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: OTU, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
samdf_out <- phyloseq::sample_data(ps2) %>%
    as("matrix") %>%
    tibble::as_tibble()

# Write out
readr::write_csv(seqtab_out, paste0("seqtab_filtered_",pcr_primers,".csv"))
readr::write_csv(taxtab_out, paste0("taxtab_filtered_",pcr_primers,".csv"))
readr::write_csv(samdf_out, paste0("samdf_filtered_",pcr_primers,".csv"))
saveRDS(ps2, paste0("ps_filtered_",pcr_primers,".rds"))



# stop(" *** stopped manually *** ") ##########################################
