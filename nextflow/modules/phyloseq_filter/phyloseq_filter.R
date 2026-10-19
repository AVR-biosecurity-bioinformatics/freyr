#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 
primers                 <- args$primers
taxtab_file             <- args$taxtab_file
seqtab_file             <- args$seqtab_file
samdf_file              <- args$samdf_file
filters_file            <- args$filters_file
fasta_file              <- args$fasta_file
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

# convert "NA" strings to true NA, or convert number strings to numeric
if(min_sample_reads %in% c(0, "NA", NA)){ min_sample_reads <- NA } else { min_sample_reads <- as.numeric(min_sample_reads) }
if(min_taxa_reads %in% c(0, "NA", NA)){ min_taxa_reads <- NA } else { min_taxa_reads <- as.numeric(min_taxa_reads) }  
if(min_taxa_ra %in% c(0, "NA", NA)){ min_taxa_ra <- NA } else { min_taxa_ra <- as.numeric(min_taxa_ra) } 

cluster_threshold <- suppressWarnings(as.integer(cluster_threshold))

### run R code

# Taxonomic filtering

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
tax_cols <- c(
    "Root",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
)

## Read input files
taxtab_dt <- data.table::fread(taxtab_file)
seqtab_dt <- data.table::fread(seqtab_file)
seq_filters <- data.table::fread(filters_file)
seqs <- Biostrings::readDNAStringSet(fasta_file)

# Restore the original ASV identifier from the annotated FASTA header.
names(seqs) <- sub("\\|.*$","",names(seqs))

# Taxtab: rows = ASVs, columns = taxonomic ranks.
taxtab <- as.matrix(taxtab_dt[, ..tax_cols])
rownames(taxtab) <- taxtab_dt$seq_name

# Seqtab: rows = ASVs, columns = samples.
sample_cols <- setdiff(names(seqtab_dt), c("seq_name", "cluster"))
seqtab <- as.matrix( seqtab_dt[, ..sample_cols] )
rownames(seqtab) <- seqtab_dt$seq_name
storage.mode(seqtab) <- "integer"

# Ensure seqtabs and taxtabs are aligned
# Require unique ASV identifiers and matching ASV sets.
stopifnot(
    !anyDuplicated(rownames(seqtab)),
    !anyDuplicated(rownames(taxtab)),
    setequal(rownames(seqtab), rownames(taxtab))
)

# Align taxonomy rows to the abundance matrix before constructing keep_taxa.
taxtab <- taxtab[rownames(seqtab), , drop = FALSE]

stopifnot(identical(rownames(seqtab), rownames(taxtab)))

# Reconstruct sample metadata with sample identifiers as row names.
samdf <- as.data.frame(data.table::fread(samdf_file))
rownames(samdf) <- samdf$sample_primers

## Per-taxon filters (removes whole ASVs)
# Convert a comma-separated taxonomic filter into a character vector.
# Examples:
#   "Insecta,Gastropoda" -> c("Insecta", "Gastropoda")
#   "NA"                 -> NULL (no filter)
#   NA_character_        -> NULL (no filter)
#   ""                   -> NULL (no filter)
parse_taxa_filter <- function(x) {
    if (is.null(x) || length(x) == 0L || all(is.na(x)) || all(trimws(x) %in% c("", "NA"))) {
        return(NULL)
    }
    taxa <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))

    # Remove empty values and literal "NA" values.
    taxa <- taxa[!is.na(taxa) & nzchar(taxa) & taxa != "NA"]

    if (length(taxa) == 0L) {
        NULL
    } else {
        unique(taxa)
    }
}

# Parse taxonomic filters. Each target may contain one or more
# comma-separated taxa, for example "Insecta,Gastropoda".
target_kingdom <- parse_taxa_filter(target_kingdom)
target_phylum  <- parse_taxa_filter(target_phylum)
target_class   <- parse_taxa_filter(target_class)
target_order   <- parse_taxa_filter(target_order)
target_family  <- parse_taxa_filter(target_family)
target_genus   <- parse_taxa_filter(target_genus)
target_species <- parse_taxa_filter(target_species)


# Return TRUE for every sequence when no filter is requested.
# Otherwise, return TRUE where the classification matches any target taxon.
matches_taxon <- function(values, targets) {
    if (is.null(targets)) {
        return(rep(TRUE, length(values)))
    }

    !is.na(values) & values %in% targets
}

# Retain only taxa matching every requested taxonomic rank.
# Multiple taxa within a rank use OR logic: Class = Insecta OR Gastropoda
# Filters across ranks use AND logic: Phylum = Arthropoda AND Class = Insecta
keep_taxa <-
    matches_taxon(taxtab[, "Kingdom"], target_kingdom) &
    matches_taxon(taxtab[, "Phylum"],  target_phylum) &
    matches_taxon(taxtab[, "Class"],   target_class) &
    matches_taxon(taxtab[, "Order"],   target_order) &
    matches_taxon(taxtab[, "Family"],  target_family) &
    matches_taxon(taxtab[, "Genus"],   target_genus) &
    matches_taxon(taxtab[, "Species"], target_species)

names(keep_taxa) <- rownames(taxtab)

# Identify sequences that passed all sequence-quality filters.
seqs_passing <- with(
    seq_filters,
    seq_name[
        chimera_filter &
        length_filter &
        phmm_filter &
        frame_filter
    ]
)

# Further restrict the retained taxa to sequences passing all soft filters.
keep_taxa <- keep_taxa & names(keep_taxa) %in% seqs_passing

# Subset the OTU matrix using retained taxa logical vector
seqtab <- seqtab[keep_taxa, , drop = FALSE]

if (nrow(seqtab) == 0L) {
    stop("No ASVs passed the taxonomic and sequence filters for primers '", primers,"'")
}

## Per-observation filters (removes observations of ASVs within samples)

# Set observations below fixed read threshold to 0
if (!is.na(min_taxa_reads)) {
    seqtab[seqtab <= min_taxa_reads] <- 0
}

# Convert observations to within-sample relative abundance, remove
# low-abundance observations, and retain the original count matrix.
if (!is.na(min_taxa_ra)) {
    sample_totals <- colSums(seqtab)
    nonzero_samples <- sample_totals > 0

    if (any(nonzero_samples)) {
        seqtab_nonzero <- seqtab[, nonzero_samples, drop = FALSE]

        ra <- sweep(seqtab_nonzero, 2, sample_totals[nonzero_samples], "/")

        seqtab_nonzero[ra <= min_taxa_ra] <- 0

        seqtab[, nonzero_samples ] <- seqtab_nonzero
    }
}

## Sample-level filters

# Remove samples below the minimum retained-read threshold.
sample_totals <- colSums(seqtab)

if (!is.na(min_sample_reads)) {
    keep_samples <- sample_totals >= min_sample_reads
    seqtab <- seqtab[, keep_samples, drop = FALSE]
} 

# Drop ASVs with zero abundance across all retained samples.
seqtab_filtered <- seqtab[rowSums(seqtab) > 0, , drop = FALSE]

if (ncol(seqtab_filtered) == 0L) {
    stop("No samples remained after filtering for primers '", primers,"'")
}

if (nrow(seqtab_filtered) == 0L) {
    stop("No ASVs remained after abundance and sample filtering for primers '", primers, "'")
}

# subset taxtab, samdf, and seqs to samples / ASVs retained through filters
taxtab_filtered <- taxtab[rownames(seqtab_filtered), ,drop = FALSE]
samdf_filtered <- samdf[colnames(seqtab_filtered), ,drop = FALSE]
seqs_filtered <- seqs[rownames(seqtab_filtered)]

# Message how many were removed
message(nrow(samdf) - nrow(samdf_filtered), " samples and ", nrow(taxtab) - nrow(taxtab_filtered), " ASVs dropped")

### cluster filtered ASVs ------------------------------------------------------------------------

if (is.na(cluster_threshold)) {

    # No clustering requested.
    cluster_dt <- data.table::data.table(
        seq_name = names(seqs_filtered),
        cluster = NA_integer_
    )

    clusters_out <- data.table::copy(cluster_dt)

} else {

    set.seed(1)

    clusters <- DECIPHER::Clusterize(
        seqs_filtered,
        cutoff = 1 - (cluster_threshold / 100),
        invertCenters = TRUE,
        # maxPhase1 = 1e5,
        # maxPhase2 = 1e4,
        # maxPhase3 = 1e4,
        # singleLinkage = FALSE,
        # rareKmers = 100,
        processors = 1
    )

    # Clusterize returns sequence names as row names and cluster assignments
    # in a one-column object.
    clusters_out <- data.table::data.table(
        seq_name = rownames(clusters),
        cluster = as.integer(clusters[, 1L])
    )

    # Internal lookup uses absolute cluster IDs.
    cluster_dt <- data.table::copy(clusters_out)
    cluster_dt[, cluster := abs(cluster)]
}

data.table::setkey(cluster_dt, seq_name)
cluster_lookup <- cluster_dt$cluster
names(cluster_lookup) <- cluster_dt$seq_name

# Write clusters_out
data.table::fwrite(clusters_out, paste0("clusters_", primers, ".csv"))

### outputs -----------------------------------------------------------------------------

# Export raw csv 
raw_filtered <- melt_phyloseq(
    otutab = seqtab_filtered,
    taxtab = taxtab_filtered,
    samdf = samdf_filtered,
    seqs = seqs_filtered
)

# Add cluster assignments
raw_filtered[, cluster := cluster_lookup[seq_name]]

# Set column order
data.table::setcolorder(
    raw_filtered,
    c(
        "seq_name",
        "sequence",
        "cluster",
        setdiff(
            names(raw_filtered),
            c(
                "seq_name",
                "sequence",
                "cluster"
            )
        )
    )
)

data.table::fwrite(raw_filtered, paste0("raw_filtered_", primers, ".csv"))

# Create summmary from phyloseq
summary_filtered <- summarise_phyloseq(
    otutab = seqtab_filtered,
    taxtab = taxtab_filtered,
    seqs = seqs_filtered
)

# Add cluster assignments
summary_filtered[, cluster := cluster_lookup[seq_name]]

# Set column order
preferred_cols <- c(
    "seq_name",
    "sequence",
    "cluster",
    tax_cols
)

data.table::setcolorder(
    summary_filtered,
    c(
        intersect(
            preferred_cols,
            names(summary_filtered)
        ),
        setdiff(
            names(summary_filtered),
            preferred_cols
        )
    )
)

data.table::fwrite(summary_filtered, paste0( "summary_filtered_",primers,".csv"))


### Phyloseq outputs -----------------------------------------------------------------------------

# save seqtab (data) as wide tibble (rows = seq_name, columns = cluster, all of sample_primers)
seqtab_out <- data.table::as.data.table(
    seqtab_filtered,
    keep.rownames = "seq_name"
)

# Add cluster assignments.
seqtab_out[, cluster := cluster_lookup[seq_name]]

# Set column order
data.table::setcolorder(seqtab_out, c( "seq_name", "cluster", setdiff( names(seqtab_out), c("seq_name", "cluster"))))

# Write out
data.table::fwrite(seqtab_out, paste0("seqtab_filtered_", primers, ".csv"))

# save taxtab as tibble (rows = seq_name, columns = taxonomic ranks (Root -> Species)) with unclassified labels converted to NA.
taxtab_out <- taxtab_filtered %>%
    unclassified_to_na() %>%
    data.table::as.data.table(
        keep.rownames = "seq_name"
    )

# Add cluster assignments.
taxtab_out[ , cluster := cluster_lookup[seq_name]]

# Place identifier columns before taxonomy columns.
data.table::setcolorder(taxtab_out, c("seq_name", "cluster", tax_cols))

# Write out
data.table::fwrite(taxtab_out, paste0("taxtab_filtered_", primers, ".csv"))

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seqs_output <- seqs_filtered

# Match taxonomy rows to the FASTA sequence order.
taxtab_fasta <- taxtab_out[ match( names(seqs_output), taxtab_out$seq_name ),]

# Collapse taxonomy ranks into a single string.
taxonomy_strings <- apply(
    taxtab_fasta[, ..tax_cols, drop = FALSE],
    1,
    paste,
    collapse = ";"
)

# Create FASTA headers.
names(seqs_output) <- paste(
    paste(
        taxtab_fasta$seq_name,
        primers,
        sep = "|"
    ),
    taxonomy_strings,
    sep = ";"
)
write_fasta(seqs_output, paste0("asvs_filtered_", primers, ".fasta"))  

# save samplesheet
data.table::as.data.table(
    samdf_filtered,
    keep.rownames = "sample_primers"
  ) %>%
    data.table::fwrite(file = paste0("samdf_filtered_",primers,".csv"))


## create phyloseq object
ps_f <- 
    phyloseq::phyloseq(
        phyloseq::otu_table(seqtab_filtered, taxa_are_rows = TRUE), 
        phyloseq::tax_table(taxtab_filtered),
        phyloseq::sample_data(samdf_filtered),
        phyloseq::refseq(seqs_filtered)
    )
saveRDS(ps_f, paste0("ps_filtered_",primers,".rds"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})