#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "phyloseq",
    "purrr",
    "readr",
    "scales",
    "stringr",
    "tibble",
    "data.table",
    "tidyr",
    "vegan",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "pcr_primers",
    "taxtab",
    "seqtab_list",
    "loci_params",
    "samdf"
)
lapply(nf_vars, nf_var_check)

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
    seqtab_final <- seqtab_list[[1]]
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
ps <- step_phyloseq(
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
# Export raw csv  - NOTE: This is memory intensive
melt_phyloseq(ps) %>%
    readr::write_csv(., paste0("raw_unfiltered_",pcr_primers,".csv"))

# Export species level summary of filtered results
summarise_phyloseq(ps) %>%
    readr::write_csv(., paste0("summary_unfiltered_",pcr_primers,".csv"))

# Output fasta of all ASVs
Biostrings::writeXStringSet(phyloseq::refseq(ps), filepath = paste0("asvs_unfiltered_",pcr_primers,".fasta"), width = 100) 

# write .nwk file if phylogeny present
if(!is.null(phyloseq::phy_tree(ps, errorIfNULL = FALSE))){
    #Output newick tree
    ape::write.tree(phyloseq::phy_tree(ps), file = paste0("tree_unfiltered",pcr_primers,".nwk"))
}

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = sample_id, cols = OTU name (hash), cells = abundance)
seqtab_out <- phyloseq::otu_table(ps) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_id")

# save taxtab as long tibble (rows = OTU/ASV, cols = tax rankings)
taxtab_out <- phyloseq::tax_table(ps) %>%
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
samdf_out <- phyloseq::sample_data(ps) %>%
    as("matrix") %>%
    tibble::as_tibble()

# Write out
readr::write_csv(seqtab_out, paste0("seqtab_unfiltered_",pcr_primers,".csv"))
readr::write_csv(taxtab_out, paste0("taxtab_unfiltered_",pcr_primers,".csv"))
readr::write_csv(samdf_out, paste0("samdf_unfiltered_",pcr_primers,".csv"))
saveRDS(ps, paste0("ps_unfiltered_",pcr_primers,".rds"))

# ## creates accumulation curve plots and saves plot
# gg.acc_curve <- rareplot(ps, step=1L, threshold = max(samdf$min_sample_reads))

# pdf(file=paste0("accumulation_curve_",pcr_primers,".pdf"), width = 11, height = 8 , paper="a4r")
#     print(gg.acc_curve)
# try(dev.off(), silent=TRUE)

# stop(" *** stopped manually *** ") ##########################################
