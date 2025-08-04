#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                 <- args$primers
ps                      <- args$ps
filters_tibble          <- args$filters_tibble
cluster_threshold       <- args$cluster_threshold
target_kingdom          <- args$target_kingdom
target_phylum           <- args$target_phylum
target_class            <- args$target_class
target_order            <- args$target_order
target_family           <- args$target_family
target_genus            <- args$target_genus
target_species          <- args$target_species
min_sample_reads        <- args$min_sample_reads
min_taxa_reads          <- args$min_taxa_reads
min_taxa_ra             <- args$min_taxa_ra

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
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables
ps <- readRDS(ps)
seq_filters <- readr::read_csv(filters_tibble)

# convert "NA" strings to true NA
if(target_kingdom == "NA"){ target_kingdom <- NA }
if(target_phylum == "NA"){ target_phylum <- NA }
if(target_class == "NA"){ target_class <- NA }
if(target_order == "NA"){ target_order <- NA }
if(target_family == "NA"){ target_family <- NA }
if(target_genus == "NA"){ target_genus <- NA }
if(target_species == "NA"){ target_species <- NA }

# convert "NA" strings to true NA, or convert number strings to numeric
if(min_sample_reads %in% c(0, "NA", NA)){ min_sample_reads <- NA } else { min_sample_reads <- as.numeric(min_sample_reads) }
if(min_taxa_reads %in% c(0, "NA", NA)){ min_taxa_reads <- NA } else { min_taxa_reads <- as.numeric(min_taxa_reads) }  
if(min_taxa_ra %in% c(0, "NA", NA)){ min_taxa_ra <- NA } else { min_taxa_ra <- as.numeric(min_taxa_ra) } 

cluster_threshold <- suppressWarnings(as.integer(cluster_threshold))

quiet <- FALSE

### run R code

## from step_filter_phyloseq()
# Taxonomic filtering
taxtab <- phyloseq::tax_table(ps) %>%
    as("matrix") %>% 
    as.data.frame()

# Check if any taxonomic filters are enabled
ps_taxfiltered <- ps

if(any(!sapply(c(target_kingdom, target_phylum, target_class, target_order, target_family, target_genus, target_species), is.na))){

    # Filter kingdom
    if(!is.na(target_kingdom)){
        if (any(stringr::str_detect(taxtab$Kingdom, target_kingdom))){
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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
        ps_taxfiltered <- ps_taxfiltered %>%
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

## remove sequences that failed the sequence-level soft filters
seqs_passing <- 
    seq_filters %>%
    dplyr::filter(
        chimera_filter, 
        length_filter,
        phmm_filter,
        frame_filter
    ) %>%
    dplyr::pull(seq_name)

ps_seqfiltered <- phyloseq::prune_taxa(taxa = seqs_passing, ps_taxfiltered)

# Remove any taxa under read count or relative abundance thresholds
### TODO: Add comments explaining what is happening here
if(!is.na(min_taxa_reads) & is.na(min_taxa_ra)){
    ps_abfiltered <- 
        phyloseq::transform_sample_counts(
            ps_seqfiltered, 
            function(OTU, ab = min_taxa_reads){ 
                ifelse(OTU <= ab,  0, OTU) 
            }
        )
} else if(is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps_abfiltered <- 
        phyloseq::transform_sample_counts(
            ps_seqfiltered, 
            function(OTU, ab = min_taxa_ra ){ 
                if (sum(OTU) == 0){
                    return(OTU)
                } else {
                    ifelse((OTU / sum(OTU)) <= ab,  0, OTU) 
                }
            }
        )
} else if (!is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps_abfiltered <- 
        ps_seqfiltered %>%
        phyloseq::transform_sample_counts(
            function(OTU, ab = min_taxa_reads){ 
                ifelse(OTU <= ab,  0, OTU) 
            }
        ) %>%
        phyloseq::transform_sample_counts(
            function(OTU, ab = min_taxa_ra ){ 
                if (sum(OTU) == 0){
                    return(OTU)
                } else {
                    ifelse((OTU / sum(OTU)) <= ab,  0, OTU) 
                }
            }
        )
} else {
    if (!quiet){message(paste0("No minimum abundance filters set - skipping this filter"))}
    ps_abfiltered <- ps_seqfiltered
}

#Remove all samples under the minimum read threshold 
if(!is.na(min_sample_reads)){

    if (all(phyloseq::sample_sums(ps_abfiltered) < min_sample_reads)) {
        stop(paste0("ERROR: No samples contained reads above the minimum threshold of ", min_sample_reads, " for primers '", primers, "' -- consider lowering this value"))
    }
    
    ps_sampfiltered <- ps_abfiltered %>%
        phyloseq::prune_samples(phyloseq::sample_sums(.) >= min_sample_reads, .) %>% 
        phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table

} else if (is.na(min_sample_reads) ) {

    if (!quiet){message(paste0("No minimum sample reads filter set - skipping this filter"))}

    ps_sampfiltered <- ps_abfiltered %>%
        phyloseq::prune_samples(phyloseq::sample_sums(.)>=0,.) %>%
        phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table

} else {
    stop(paste0("--min_sample_reads value is '",min_sample_reads,"', which is not allowed."))
}

# Message how many were removed
if(!quiet){message(phyloseq::nsamples(ps) - phyloseq::nsamples(ps_sampfiltered), " Samples and ", phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps_sampfiltered), " ASVs dropped")}

### cluster filtered sequences

seqs <- ps_sampfiltered %>% phyloseq::refseq()

if (cluster_threshold %>% is.na){
  
    # make clusters "NA" if threshold not given   
    cluster_tibble <- 
        tibble::as_tibble_col(names(seqs_combined), column_name = "seq_name") %>%
        dplyr::mutate(cluster = NA_integer_)

    # for export
    clusters_out <- cluster_tibble

} else {
    set.seed <- 1; clusters <- DECIPHER::Clusterize(
        seqs, 
        cutoff = 1 - (cluster_threshold / 100),
        invertCenters = TRUE,
        # maxPhase1 = 1e5,
        # maxPhase2 = 1e4,
        # maxPhase3 = 1e4, 
        # singleLinkage = FALSE,
        # rareKmers = 100,
        processors = 1
    )

    # for internal use
    cluster_tibble <- 
        clusters %>% 
        tibble::as_tibble(rownames = "seq_name") %>%
        dplyr::mutate(cluster = abs(cluster))

    # for export
    clusters_out <- 
        clusters %>% 
        tibble::as_tibble(rownames = "seq_name") 
  
}

### output filtered results per primer pair; from step_output_summary()

# Export raw csv  - NOTE: This is memory intensive
melt_phyloseq(ps_sampfiltered) %>% 
    tibble::as_tibble() %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = sequence) %>%
    readr::write_csv(., paste0("raw_filtered_",primers,".csv"))

# Export species level summary of filtered results
summarise_phyloseq(ps_sampfiltered) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>%
    dplyr::relocate(cluster, .after = sequence) %>%
    readr::write_csv(., paste0("summary_filtered_",primers,".csv"))

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seqs_output <- phyloseq::refseq(ps_sampfiltered)

seq_names_new <- 
    phyloseq::tax_table(ps_sampfiltered) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # ensure sequence name order is same as the DSS object
    dplyr::arrange(factor(seq_name, levels = names(seqs_output))) %>%
    # add primers
    dplyr::mutate(primers = primers, .after = seq_name) %>%
    # unite columns into a single header string per sequence
    tidyr::unite(col = "lineage", Root:Species, sep = ";") %>%
    tidyr::unite(col = "id", c(seq_name, primers), sep = "|") %>%
    tidyr::unite(col = "header", c(id, lineage), sep = ";") %>%
    dplyr::pull(header)

names(seqs_output) <- seq_names_new

write_fasta(seqs_output, paste0("asvs_unfiltered_", primers, ".fasta"))  

## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = seq_name, cols = ASV name (hash), cells = abundance)
seqtab_out <- phyloseq::otu_table(ps_sampfiltered) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = seq_name)

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out <- phyloseq::tax_table(ps_sampfiltered) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = seq_name) %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    )

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
if(!all(colnames(taxtab_out) == c("seq_name", "cluster", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
            Database requires the columns: seq_name, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
}

# save samplesheet
samdf_out <- phyloseq::sample_data(ps_sampfiltered) %>%
    as("matrix") %>%
    tibble::as_tibble()

# Write out
readr::write_csv(clusters_out, paste0("clusters_",primers,".csv"))
readr::write_csv(seqtab_out, paste0("seqtab_filtered_",primers,".csv"))
readr::write_csv(taxtab_out, paste0("taxtab_filtered_",primers,".csv"))
readr::write_csv(samdf_out, paste0("samdf_filtered_",primers,".csv"))
saveRDS(ps_sampfiltered, paste0("ps_filtered_",primers,".rds"))


# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})