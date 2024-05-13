#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    # "bs4Dash",
    # "clustermq",
    "dada2",
    # "DECIPHER",
    "dplyr",
    # "future",
    # "ggplot2",
    # "gridExtra",
    # "gt",
    # "magrittr",
    # "markdown",
    # "ngsReports",
    # "patchwork",
    "phyloseq",
    # "pingr",
    # "purrr",
    "readr",
    # "rlang",
    # "rstudioapi",
    # "savR",
    # "scales",
    # "seqateurs",
    # "shiny",
    # "shinybusy",
    # "shinyWidgets",
    # "ShortRead",
    "stringr",
    # "taxreturn",
    "tibble",
    # "tidyr",
    # "vegan",
    # "visNetwork",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))


## check and define variables
taxtab <- readRDS(taxtab)
taxtab %>% # save for debugging
    tibble::as_tibble(rownames = "OTU") %>%
    readr::write_csv(., paste0("taxtab_", pcr_primers, ".csv")) 

seqtab_list <- # convert Groovy to R list format
    stringr::str_extract_all(seqtab_list, pattern = "[^\\s,\\[\\]]+") %>% unlist()
seqtab_list <- lapply(seqtab_list, readRDS) # read in seqtabs and store as list of matrices

samdf <- readr::read_csv(samdf, show_col_types = FALSE)

### run R code

## merge sequence tables across flowcells and loci
if ( length(seqtab_list) > 1 ){ # if there is more than one seqtab, merge together
    seqtab_final <- dada2::mergeSequenceTables(tables=seqtab_list)
} else if( length(seqtab_list) == 1 ) { # if there is only one seqtab, keep and unlist
    seqtab_final <- seqtab_list %>% unlist()
}
seqtab_final %>% # save for debugging
    tibble::as_tibble(rownames = "OTU") %>% 
    readr::write_csv(., paste0("seqtab_final_", pcr_primers, ".csv"))

## mutate samdf to add pcr_primers to sample_id, to make consistent with new seqtab format
samdf_renamed <- samdf %>% 
    dplyr::mutate(
        sample_id_orig = sample_id, # save original sample_id as a new column
        sample_id = paste0(sample_id,"_",pcr_primers) # mutate sample_id to add primer id at the end
        )

## run step_phyloseq() on merged seqtabs, taxtabs and samplesheet to create phyloseq object
ps <-   step_phyloseq(
            seqtab = seqtab_final,
            taxtab = taxtab,
            samdf = samdf_renamed,
            seqs = NULL,
            phylo = NULL,
            name_variants = FALSE
        )

## name OTUs using hash
phyloseq::taxa_names(ps) <- phyloseq::tax_table(ps)[,ncol(phyloseq::tax_table(ps))] # use final column of tax_table (hash) to name OTUs
phyloseq::tax_table(ps) <- phyloseq::tax_table(ps)[,1:ncol(phyloseq::tax_table(ps))-1] # remove hash 'rank' from tax_table

## output summaries; from step_output_summary()
# Export raw csv
phyloseq::psmelt(ps) %>% 
    filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    write_csv(., paste0("raw_unfiltered_",pcr_primers,".csv"))

# Export species level summary of filtered results
phyloseq::psmelt(ps) %>% 
    filter(Abundance > 0) %>%
    left_join(
        refseq(ps) %>% as.character() %>% enframe(name="OTU", value="sequence"),
        by = "OTU"
        ) %>%
    dplyr::select(OTU, sequence, rank_names(ps), sample_id, Abundance ) %>%
    pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    write_csv(., paste0("summary_unfiltered_",pcr_primers,".csv"))

# Output fasta of all ASVs
seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps)))
Biostrings::writeXStringSet(seqs, filepath = paste0("asvs_unfiltered_",pcr_primers,".fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phy_tree(ps, errorIfNULL = FALSE))){
    #Output newick tree
    write.tree(phy_tree(ps), file = paste0("tree_unfiltered",pcr_primers,".nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out <- phyloseq::otu_table(ps) %>%
    as("matrix") %>%
    as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out <- phyloseq::tax_table(ps) %>%
    as("matrix") %>%
    as_tibble(rownames = "OTU") %>%
    seqateurs::unclassified_to_na(rownames = FALSE)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab) == c("OTU", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: OTU, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
samdf_out <- phyloseq::sample_data(ps) %>%
    as("matrix") %>%
    as_tibble()

# Write out
write_csv(seqtab_out, paste0("seqtab_unfiltered_",pcr_primers,".csv"))
write_csv(taxtab_out, paste0("taxtab_unfiltered_",pcr_primers,".csv"))
write_csv(samdf_out, paste0("samdf_unfiltered_",pcr_primers,".csv"))
saveRDS(ps, paste0("ps_unfiltered_",pcr_primers,".rds"))

## creates accumulation curve plots and saves plot
gg.acc_curve <- rareplot(ps, step="auto", threshold = max(samdf$min_sample_reads))

pdf(file=paste0("accumulation_curve_",pcr_primers,".pdf"), width = 11, height = 8 , paper="a4r")
    print(gg.acc_curve)
try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################
