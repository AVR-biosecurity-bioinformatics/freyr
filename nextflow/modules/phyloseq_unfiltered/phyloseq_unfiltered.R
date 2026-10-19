#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 
primers                     <- args$primers
taxtab_file                 <- args$taxtab_file
seqtab_file                 <- args$seqtab_file
filters_file                <- args$filters_file
fasta_file                  <- args$fasta_file
samplesheet_split_file      <- args$samplesheet_split_file
sample_metadata_file        <- args$sample_metadata_file
cluster_threshold           <- args$cluster_threshold

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
    "ggplot2",
    "phyloseq",
    "readr",
    "scales",
    "stringr",
    "tibble",
    "tidyr",
    "data.table",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

## Read in objects
taxtab <- data.table::fread(taxtab_file)
seqtab_combined <- data.table::fread(seqtab_file)
filters <- data.table::fread(filters_file)
seqs <- Biostrings::readDNAStringSet(fasta_file)
cluster_threshold <- suppressWarnings(as.integer(cluster_threshold))

# Create the complete sample metadata table.
samples_df <- data.table::fread(samplesheet_split_file)[
  data.table::fread(sample_metadata_file),
  on = .(sample, read_group)
] %>%
  as.data.frame() %>%
  dplyr::filter(primers == .env$primers)

rownames(samples_df) <- samples_df$sample_primers

# create phyloseq-format otu_table (seqtab), ignoring filters -- ASV names are rows
otu_cols <- setdiff( names(seqtab_combined), c("seq_name", "sequence", "cluster"))
otutab_matrix <- as.matrix(seqtab_combined[, ..otu_cols])
rownames(otutab_matrix) <- seqtab_combined[["seq_name"]]
storage.mode(otutab_matrix) <- "integer"

# create phyloseq-format taxtab
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
taxtab_matrix <- as.matrix(taxtab[, ..tax_cols])
rownames(taxtab_matrix) <- taxtab[["seq_name"]]

# validate all objects are compatible with phyloseq
validate_phyloseq_inputs(
    otutab = otutab_matrix,
    taxtab = taxtab_matrix,
    samdf = samples_df,
    seqs = seqs,
    tax_cols = tax_cols
)

### Cluster unfiltered ASVs -----------------------------------------------

if (is.na(cluster_threshold)) {

    # No clustering requested.
    cluster_dt <- data.table::data.table(
        seq_name = names(seqs),
        cluster = NA_integer_
    )

    clusters_out <- data.table::copy(cluster_dt)

} else {

    set.seed(1)

    clusters <- DECIPHER::Clusterize(
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

# Export raw CSV. No longer requires phyloseq object
raw_dt <- melt_phyloseq(
    otutab = otutab_matrix,
    taxtab = taxtab_matrix,
    samdf = samples_df,
    seqs = seqs
)

# Add cluster assignments.
raw_dt[, cluster := cluster_lookup[seq_name]]

# Set column order
data.table::setcolorder(
    raw_dt,
    c(
        "seq_name",
        "sequence",
        "cluster",
        setdiff(
            names(raw_dt),
            c("seq_name", "sequence", "cluster")
        )
    )
)

data.table::fwrite(
    raw_dt,
    paste0(
        "raw_unfiltered_",
        primers,
        ".csv"
    )
)


# Export summary CSV. No longer requires phyloseq object
summary_dt <- summarise_phyloseq(
    otutab = otutab_matrix,
    taxtab = taxtab_matrix,
    seqs = seqs
)

summary_dt <- filters[
    summary_dt,
    on = c(
        "seq_name",
        "sequence"
    )
]

# Add cluster assignments.
summary_dt[, cluster := cluster_lookup[seq_name]]

preferred_cols <- c(
    "seq_name",
    "sequence",
    "cluster",
    tax_cols,
    "chimera_filter",
    "length_filter",
    "phmm_filter",
    "frame_filter"
)

preferred_cols <- intersect(
    preferred_cols,
    names(summary_dt)
)

# Set column order
data.table::setcolorder(
    summary_dt,
    c(
        preferred_cols,
        setdiff(
            names(summary_dt),
            preferred_cols
        )
    )
)

data.table::fwrite(summary_dt, paste0("summary_unfiltered_", primers, ".csv"))

# save seqtab (data) as wide tibble (rows = seq_name, columns = cluster, all of sample_primers)
seqtab_out <- data.table::as.data.table(
    otutab_matrix,
    keep.rownames = "seq_name"
)

# Add cluster assignments.
seqtab_out[, cluster := cluster_lookup[seq_name]]

# Set column order
data.table::setcolorder(seqtab_out, c( "seq_name", "cluster", setdiff( names(seqtab_out), c("seq_name", "cluster"))))

# Write out
data.table::fwrite(seqtab_out, paste0("seqtab_unfiltered_", primers, ".csv"))

# save taxtab as tibble (rows = seq_name, columns = taxonomic ranks (Root -> Species)) with unclassified labels converted to NA.
# Use taxtab_matrix directy, no need to re-extract from phyloseq
taxtab_out <- taxtab_matrix %>%
    unclassified_to_na() %>%
    data.table::as.data.table(
        keep.rownames = "seq_name"
    )

# Add cluster assignments.
taxtab_out[ , cluster := cluster_lookup[seq_name]]

# Place identifier columns before taxonomy columns.
data.table::setcolorder(taxtab_out, c("seq_name", "cluster", tax_cols))

# Write out
data.table::fwrite(taxtab_out, paste0("taxtab_unfiltered_", primers, ".csv"))

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seqs_output <- seqs

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
write_fasta(seqs_output, paste0("asvs_unfiltered_", primers, ".fasta"))  

# save samplesheet
data.table::as.data.table(
    samples_df,
    keep.rownames = "sample_primers"
  ) %>%
    data.table::fwrite(file = paste0("samdf_unfiltered_",primers,".csv"))


## create phyloseq object
ps_uf <- 
    phyloseq::phyloseq(
        phyloseq::otu_table(otutab_matrix, taxa_are_rows = TRUE), 
        phyloseq::tax_table(taxtab_matrix),
        phyloseq::sample_data(samples_df),
        phyloseq::refseq(seqs)
    )
saveRDS(ps_uf, paste0("ps_unfiltered_",primers,".rds"))

# Write filters out
data.table::fwrite(filters,paste0("filters_", primers, ".csv"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})