#!/usr/bin/env Rscript
.libPaths(.Library)

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

### run R code

# Taxonomic filtering

# Extract the taxonomy and OTU table once.
taxtab <- as(phyloseq::tax_table(ps), "matrix")
otutab <- as(phyloseq::otu_table(ps), "matrix")
taxa_are_rows <- phyloseq::taxa_are_rows(ps)

if (taxa_are_rows) {
    otutab <- t(otutab)
}

# Logical vector of taxa
keep_taxa <- rep(TRUE, nrow(taxtab))
names(keep_taxa) <- rownames(taxtab)

# Taxonomic filters
if (!is.na(target_kingdom)) keep_taxa <- keep_taxa & taxtab[, "Kingdom"] == target_kingdom
if (!is.na(target_phylum))  keep_taxa <- keep_taxa & taxtab[, "Phylum"]  == target_phylum
if (!is.na(target_class))   keep_taxa <- keep_taxa & taxtab[, "Class"]   == target_class
if (!is.na(target_order))   keep_taxa <- keep_taxa & taxtab[, "Order"]   == target_order
if (!is.na(target_family))  keep_taxa <- keep_taxa & taxtab[, "Family"]  == target_family
if (!is.na(target_genus))   keep_taxa <- keep_taxa & taxtab[, "Genus"]   == target_genus
if (!is.na(target_species)) keep_taxa <- keep_taxa & taxtab[, "Species"] == target_species

# Sequence level soft filters
seqs_passing <- with(
    seq_filters,
    seq_name[
        chimera_filter &
        length_filter &
        phmm_filter &
        frame_filter
    ]
)

keep_taxa <- keep_taxa & names(keep_taxa) %in% seqs_passing
otutab <- otutab[, keep_taxa, drop = FALSE]

if (ncol(otutab) == 0) {
    stop("No ASVs passed the taxonomic and sequence filters for primers '",primers,"'")
}

# Fixed read threshold
if (!is.na(min_taxa_reads)) {
    otutab[otutab <= min_taxa_reads] <- 0
}

# Relative abundance threshold
if (!is.na(min_taxa_ra)) {
    sample_totals <- rowSums(otutab)
    nonzero <- sample_totals > 0

    if (any(nonzero)) {
        otutab_nonzero <- otutab[nonzero, , drop = FALSE]

        ra <- sweep(
            otutab_nonzero,
            1,
            sample_totals[nonzero],
            "/"
        )

        otutab_nonzero[ra <= min_taxa_ra] <- 0
        otutab[nonzero, ] <- otutab_nonzero

    }
}

# Remove samples below the minimum retained-read threshold.
sample_totals <- rowSums(otutab)

if (!is.na(min_sample_reads)) {
    keep_samples <- sample_totals >= min_sample_reads

    if (!any(keep_samples)) {
        stop(
            "No samples contained reads above the minimum threshold of ",
            min_sample_reads,
            " for primers '", primers,
            "' -- consider lowering this value"
        )
    }
    otutab <- otutab[keep_samples, , drop = FALSE]
} 

# Drop taxa with zero abundance across all retained samples
otutab <- otutab[, colSums(otutab) > 0, drop = FALSE]

if (ncol(otutab) == 0) {
    stop("No ASVs remained after abundance and sample filtering for primers '",primers,"'")
}

# Prune the phyloseq object once for samples and once for taxa so that all
# components, including taxonomy and reference sequences, remain synchronized.
ps_sampfiltered <- phyloseq::prune_taxa(
    colnames(otutab),
    ps
)

ps_sampfiltered <- phyloseq::prune_samples(
    rownames(otutab),
    ps_sampfiltered
)

# Replace the abundance table with the already-filtered matrix.
if (taxa_are_rows ) {
    phyloseq::otu_table(ps_sampfiltered) <- phyloseq::otu_table(
        t(otutab),
        taxa_are_rows = TRUE
    )
} else {
    phyloseq::otu_table(ps_sampfiltered) <- phyloseq::otu_table(
        otutab,
        taxa_are_rows = FALSE
    )
}

# Message how many were removed
message(phyloseq::nsamples(ps) - phyloseq::nsamples(ps_sampfiltered), " Samples and ", phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps_sampfiltered), " ASVs dropped")


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
    set.seed(1); clusters <- DECIPHER::Clusterize(
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
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = sequence) %>%
    readr::write_csv(., paste0("raw_filtered_",primers,".csv"))

# Export species level summary of filtered results
summarise_phyloseq(ps_sampfiltered) %>%
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

write_fasta(seqs_output, paste0("asvs_filtered_", primers, ".fasta"))  

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