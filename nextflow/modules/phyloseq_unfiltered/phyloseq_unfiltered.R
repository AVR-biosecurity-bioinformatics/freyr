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

## check and define variables
taxtab <- data.table::fread(taxtab_file)

seqtab_combined <- data.table::fread(seqtab_file)

filters <- data.table::fread(filters_file)

seqs <- Biostrings::readDNAStringSet(fasta_file)

samplesheet_split <- readr::read_csv(samplesheet_split_file, show_col_types = FALSE)

sample_metadata <- readr::read_csv(sample_metadata_file, show_col_types = FALSE)

cluster_threshold <- suppressWarnings(as.integer(cluster_threshold))

## check taxonomy, seqtab, filters, and fasta all have the same sequences, none are missing
if(!setequal(seqtab_combined$seq_name, taxtab$seq_name)){
    stop("seqtab_combined and taxtab do not contain the exact same sequence names")
}

if(!setequal(seqtab_combined$seq_name, filters$seq_name)){
    stop("seqtab_combined and filters do not contain the exact same sequence names")
} 

if(!setequal(seqtab_combined$seq_name, names(seqs))){
    stop("seqtab_combined and seqs do not contain the exact same sequence names")
}

# Validate inported taxonomy 

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

observed_tax_cols <- setdiff(names(taxtab), "seq_name")

if (!identical(observed_tax_cols, tax_cols)) {
  warning(
    "Taxonomy ranks do not meet staging database expectations.\n",
    "Expected: ", paste(tax_cols, collapse = ", "), "\n",
    "Observed: ", paste(observed_tax_cols, collapse = ", ")
  )
}

### run R code

## create complete samplesheet with arbitrary metadata from original --samplesheet .csv
samplesheet_complete <- 
    samplesheet_split %>%
    dplyr::left_join(., sample_metadata, by = dplyr::join_by("sample","read_group"))

# create phyloseq-format otu_table (seqtab), ignoring filters -- ASV names are rows
otu_cols <- setdiff(
  names(seqtab_combined),
  c("seq_name", "sequence", "cluster")
)

otutab_matrix <- as.matrix(
  seqtab_combined[, ..otu_cols]
)

rownames(otutab_matrix) <- seqtab_combined[["seq_name"]]

missing_samples <- setdiff(
  colnames(otutab_matrix),
  samplesheet_complete[["sample_primers"]]
)

if (length(missing_samples)) {
  stop(
    "OTU-table samples are missing from the sample metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Force integer storage for discrete counts.
storage.mode(otutab_matrix) <- "integer"

otutab_ps <- phyloseq::otu_table(
  otutab_matrix,
  taxa_are_rows = TRUE
)

# create dada2-format taxtab
taxtab_matrix <- as.matrix(
  taxtab[, ..tax_cols]
)

rownames(taxtab_matrix) <- taxtab[["seq_name"]]

taxtab_ps <- phyloseq::tax_table(
  taxtab_matrix
)

# create sample information phyloseq object
samples_df <- as.data.frame(
  samplesheet_complete
)

if (anyDuplicated(samples_df[["sample_primers"]])) {
  stop("samplesheet_complete contains duplicated sample_primers values.")
}

rownames(samples_df) <- samples_df[["sample_primers"]]

samples_ps <- phyloseq::sample_data(
  samples_df
)

# create refseq object
refseq_ps <- 
    seqs %>%
    phyloseq::refseq()


## create phyloseq object
ps_uf <- 
    phyloseq::phyloseq(
        otutab_ps, 
        taxtab_ps,
        samples_ps,
        refseq_ps
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

### outputs -----------------------------------------------------------------------------

# Export raw csv 
write_raw_phyloseq(
    ps = ps_uf,
    cluster_dt = cluster_dt,
    filename = paste0(
        "raw_unfiltered_",
        primers,
        ".csv"
    )
)

# Export summary CSV
write_phyloseq_summary(
    ps = ps_uf,
    cluster_dt = cluster_dt,
    filters = filters,
    tax_cols = tax_cols,
    filename = paste0(
        "summary_unfiltered_",
        primers,
        ".csv"
    )
)

### Phyloseq outputs -----------------------------------------------------------------------------

# save seqtab (data) as wide tibble (rows = seq_name, columns = cluster, all of sample_primers)
seqtab_out <- data.table::as.data.table(
  as(phyloseq::otu_table(ps_uf), "matrix"),
  keep.rownames = "seq_name"
)

seqtab_out <- cluster_dt[
  seqtab_out,
  on = "seq_name"
]

data.table::setcolorder(
  seqtab_out,
  c(
    "seq_name",
    "cluster",
    setdiff(
      names(seqtab_out),
      c("seq_name", "cluster")
    )
  )
)

# save taxtab as tibble (rows = seq_name, columns = taxonomic ranks (Root -> Species))
taxtab_out <- data.table::as.data.table(
  as(phyloseq::tax_table(ps_uf), "matrix"),
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
    tax_cols
  )
)

# Convert propagated taxonomy labels to NA.
for (column in tax_cols) {
  data.table::set(
    taxtab_out,
    i = which(
      !is.na(taxtab_out[[column]]) &
        grepl("^\\w__", taxtab_out[[column]])
    ),
    j = column,
    value = NA_character_
  )
}

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seqs_output <- seqs

taxtab_fasta <- taxtab_out[
  match(names(seqs_output), seq_name)
]

if (anyNA(taxtab_fasta[["seq_name"]])) {
  stop("Could not match all FASTA sequences to the taxonomy table.")
}

seq_names_new <- taxtab_fasta[
  ,
  paste(
    paste(seq_name, primers, sep = "|"),
    do.call(
      paste,
      c(.SD, sep = ";")
    ),
    sep = ";"
  ),
  .SDcols = tax_cols
]

names(seqs_output) <- seq_names_new

# save samdf 
samdf_out <- data.table::as.data.table(
  as(phyloseq::sample_data(ps_uf), "data.frame"),
  keep.rownames = "sample_id"
)

if ("sample_primers" %in% names(samdf_out)) {
  if (!identical(
    as.character(samdf_out[["sample_id"]]),
    as.character(samdf_out[["sample_primers"]])
  )) {
    stop(
      "sample_data row names do not match its sample_primers column."
    )
  }

  samdf_out[, sample_id := NULL]
} else {
  data.table::setnames(
    samdf_out,
    "sample_id",
    "sample_primers"
  )
}

# Write out
data.table::fwrite(clusters_out, paste0("clusters_", primers, ".csv"))
data.table::fwrite(filters,paste0("filters_", primers, ".csv"))
data.table::fwrite(seqtab_out, paste0("seqtab_unfiltered_", primers, ".csv"))
data.table::fwrite(taxtab_out, paste0("taxtab_unfiltered_", primers, ".csv"))
data.table::fwrite(samdf_out, paste0("samdf_unfiltered_", primers, ".csv"))
write_fasta(seqs_output, paste0("asvs_unfiltered_", primers, ".fasta"))  
saveRDS(ps_uf, paste0("ps_unfiltered_",primers,".rds"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})