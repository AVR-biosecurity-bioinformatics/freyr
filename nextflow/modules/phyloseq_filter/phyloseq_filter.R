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
otutab <- otutab[, keep_taxa, drop = FALSE]

if (ncol(otutab) == 0) {
    stop("No ASVs passed the taxonomic and sequence filters for primers '",primers,"'")
}

## Per-observation filters (removes observations of ASVs within samples)

# Set observations below fixed read threshold to 0
if (!is.na(min_taxa_reads)) {
    otutab[otutab <= min_taxa_reads] <- 0
}

# Convert to relative abundance, remove low-abundance observations,
# and retain the count matrix.
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

## Sample-level filters

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
cluster_lookup <- cluster_dt$cluster
names(cluster_lookup) <- cluster_dt$seq_name

### outputs -----------------------------------------------------------------------------

# Export raw csv 
# TODO: melt_phyloseq now handles raw seqtab, taxtab, samdf - update to use these
raw_filtered <- cluster_dt[
    melt_phyloseq(
        ps = ps_sampfiltered
    ),
    on = "seq_name"
]

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
# TODO: melt_phyloseq now handles raw seqtab, taxtab, samdf - update to use these
summary_filtered <- summarise_phyloseq(ps_sampfiltered)

summary_filtered[
    ,
    cluster := cluster_lookup[seq_name]
]

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


## output phyloseq and component data; from step_output_ps

# save seqtab as wide tibble (rows = seq_name, cols = ASV name (hash), cells = abundance)
seqtab_out <- data.table::as.data.table(
    as(phyloseq::otu_table(ps_sampfiltered), "matrix"),
    keep.rownames = "seq_name"
)
seqtab_out <- cluster_dt[seqtab_out, on = "seq_name"]
data.table::setcolorder(seqtab_out, c("seq_name", "cluster", setdiff(names(seqtab_out), c("seq_name", "cluster"))) )

# save taxtab as long tibble (rows = ASV, cols = tax rankings) with unclassified labels converted to NA.
taxtab_out <- phyloseq::tax_table(ps_sampfiltered) %>%
  as("matrix") %>%
  unclassified_to_na() %>%
  data.table::as.data.table(
    keep.rownames = "seq_name"
  )

# Add cluster assignments.
taxtab_out <- cluster_dt[
  taxtab_out,
  on = "seq_name"
]

# Place identifier columns before the taxonomy columns.
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

# Match taxonomy rows to the FASTA sequence order.
taxtab_fasta <- taxtab_out[
    match(
        names(seqs_output),
        taxtab_out$seq_name
    ),
]

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