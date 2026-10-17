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

if (!identical(phyloseq::rank_names(ps), tax_cols)) {
    warning(
        "Taxonomy ranks do not meet staging database expectations.\n",
        "Expected: ", paste(tax_cols, collapse = ", "), "\n",
        "Observed: ", paste(phyloseq::rank_names(ps), collapse = ", ")
    )
}

# Extract the taxonomy and OTU table once.
taxtab <- as(phyloseq::tax_table(ps), "matrix")
otutab <- as(phyloseq::otu_table(ps), "matrix")
storage.mode(otutab) <- "integer" # Force integer storage of discrete count matrix
taxa_are_rows <- phyloseq::taxa_are_rows(ps)

if (taxa_are_rows) {
    otutab <- t(otutab)
}

# Logical vector of taxa
keep_taxa <- rep(TRUE, nrow(taxtab))
names(keep_taxa) <- rownames(taxtab)

# Helper function that handles missing taxonomy
matches_taxon <- function(values, target) {
    is.na(target) | (!is.na(values) & values == target)
}

# Taxonomic filters
keep_taxa <- keep_taxa &
    matches_taxon(taxtab[, "Kingdom"], target_kingdom) &
    matches_taxon(taxtab[, "Phylum"], target_phylum) &
    matches_taxon(taxtab[, "Class"], target_class) &
    matches_taxon(taxtab[, "Order"], target_order) &
    matches_taxon(taxtab[, "Family"], target_family) &
    matches_taxon(taxtab[, "Genus"], target_genus) &
    matches_taxon(taxtab[, "Species"], target_species)

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


### cluster filtered ASVs ------------------------------------------------------------------------

seqs <- ps_sampfiltered %>% phyloseq::refseq()

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

### outputs -----------------------------------------------------------------------------

# Export raw csv 
write_raw_phyloseq(
    ps = ps_sampfiltered,
    cluster_dt = cluster_dt,
    filename = paste0(
        "raw_filtered_",
        primers,
        ".csv"
    )
)

# Create summmary from phyloseq
write_phyloseq_summary(
    ps = ps_sampfiltered,
    cluster_dt = cluster_dt,
    tax_cols = tax_cols,
    filename = paste0(
        "summary_filtered_",
        primers,
        ".csv"
    )
)

### Phyloseq outputs -----------------------------------------------------------------------------


## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = seq_name, cols = ASV name (hash), cells = abundance)
seqtab_out <- data.table::as.data.table(
    as(phyloseq::otu_table(ps_sampfiltered), "matrix"),
    keep.rownames = "seq_name"
)
seqtab_out <- cluster_dt[seqtab_out, on = "seq_name"]
data.table::setcolorder(seqtab_out, c("seq_name", "cluster", setdiff(names(seqtab_out), c("seq_name", "cluster"))) )

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out <- data.table::as.data.table(
    as(phyloseq::tax_table(ps_sampfiltered), "matrix"),
    keep.rownames = "seq_name"
)

taxtab_out <- cluster_dt[
    taxtab_out,
    on = "seq_name"
]

data.table::setcolorder(
    taxtab_out,
    c(
        "seq_name",
        "cluster",
        setdiff(
            names(taxtab_out),
            c("seq_name", "cluster")
        )
    )
)

for (col in tax_cols) {
    taxtab_out[
        grepl("^\\w__", get(col)),
        (col) := NA_character_
    ]
}

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
expected_cols <- c("seq_name", "cluster", tax_cols)
if (!identical(names(taxtab_out), expected_cols)) {
    message(
        "Warning: Taxonomy table columns do not meet expectations for the staging database.\n",
        "Database requires the columns: ", paste(expected_cols, collapse = ", ")
    )
}

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seqs_output <- phyloseq::refseq(ps_sampfiltered)

# Reuse taxtab output
taxtab_fasta <- taxtab_out[
    match(names(seqs_output), seq_name)
]

seq_names_new <- taxtab_fasta[
    ,
    paste(
        paste(seq_name, primers, sep = "|"),
        do.call(paste, c(.SD, sep = ";")),
        sep = ";"
    ),
    .SDcols = tax_cols
]

names(seqs_output) <- seq_names_new

# save samplesheet
samdf_out <- data.table::as.data.table(
    as(phyloseq::sample_data(ps_sampfiltered), "data.frame"),
    keep.rownames = "sample_primers"
)

# Write out
data.table::fwrite(clusters_out, paste0("clusters_", primers, ".csv"))
data.table::fwrite(seqtab_out, paste0("seqtab_filtered_", primers, ".csv"))
data.table::fwrite(taxtab_out, paste0("taxtab_filtered_", primers, ".csv"))
data.table::fwrite(samdf_out, paste0("samdf_filtered_", primers, ".csv"))
write_fasta(seqs_output, paste0("asvs_filtered_", primers, ".fasta"))  
saveRDS(ps_sampfiltered, paste0("ps_filtered_",primers,".rds"))


# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})