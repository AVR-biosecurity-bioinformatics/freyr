#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

ps_unfiltered               <- args$ps_unfiltered
ps_filtered                 <- args$ps_filtered
unfiltered_fastas           <- args$unfiltered_fastas
filtered_fastas             <- args$filtered_fastas

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


### run R code

### unfiltered ---------------------------------------------------------------------------------------------

## merge unfiltered phyloseq objects - These must have disjoint samples
ps_u <- ps_unfiltered %>%
  stringr::str_extract_all(pattern = "[^\\s,\\[\\]]+") %>%
  unlist() %>%
  lapply(readRDS) %>%
  merge_phyloseq_disjoint_samples()

# Merge DNAStringSet objects, retaining the first occurrence of each ASV name, and write out
unfiltered_fastas %>%
  stringr::str_split_1(pattern = " ") %>%
  lapply(Biostrings::readDNAStringSet) %>%
  merge_refseq() %>%
  write_fasta("asvs_unfiltered.fasta")

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

if (!identical(phyloseq::rank_names(ps_u), tax_cols)) {
    warning(
        "Taxonomy ranks do not meet staging database expectations.\n",
        "Expected: ", paste(tax_cols, collapse = ", "), "\n",
        "Observed: ", paste(phyloseq::rank_names(ps_u), collapse = ", ")
    )
}

## output merged data tables 
data.table::fwrite(
    melt_phyloseq(ps_u),
    "raw_unfiltered.csv"
)

# Export species level summary of unfiltered results
ps_u_summary <- summarise_phyloseq(ps_u)
data.table::fwrite(ps_u_summary, "summary_unfiltered.csv")

## output phyloseq and component data

# Write out whole phyloseq object
saveRDS(ps_u, paste0("ps_unfiltered.rds"))

## read tracking output from unfiltered phyloseq
## Summarise classified reads at each taxonomic rank.
samdf <- data.table::as.data.table(
  as(phyloseq::sample_data(ps_u), "data.frame"),
  keep.rownames = "sample_primers"
)

sample_cols <- samdf$sample_primers

abundance <- as.matrix(ps_u_summary[, ..sample_cols])
storage.mode(abundance) <- "integer"

readsout <- data.table::rbindlist(
  lapply(tax_cols, function(rank) {
    data.table::data.table(
      stage = paste0("classified_", tolower(rank)),
      sample_primers = sample_cols,
      pairs = colSums(
        abundance[
          !is.na(ps_u_summary[[rank]]) &
            !grepl("^\\w__", ps_u_summary[[rank]]),
          ,
          drop = FALSE
        ]
      )
    )
  })
)

readsout <- samdf[
  readsout,
  on = "sample_primers"
][
  ,
  .(stage, sample_primers, read_group, primers, pairs)
]

data.table::fwrite(readsout,"ps_u_readsout.csv" )

# Clean up
rm(abundance, readsout, ps_u_summary)
gc()

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    data.table::as.data.table(
        keep.rownames = "sample_primers"
    ) %>%
    data.table::fwrite(file = "seqtab_unfiltered.csv")

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    unclassified_to_na() %>%
    data.table::as.data.table(keep.rownames = "seq_name") %>%
    data.table::fwrite(file = "taxtab_unfiltered.csv")

# save samplesheet
phyloseq::sample_data(ps_u) %>%
    as("data.frame") %>%
    data.table::as.data.table(keep.rownames = "sample_primers") %>%
    data.table::fwrite(file = "samdf_unfiltered.csv")

# Clean up
rm(ps_u)
gc()

### filtered ----------------------------------------------------------------------------------------------------------

## Merge filtered phyloseq objects with disjoint samples.
ps_f <- ps_filtered %>%
  stringr::str_extract_all(pattern = "[^\\s,\\[\\]]+") %>%
  unlist() %>%
  lapply(readRDS) %>%
  merge_phyloseq_disjoint_samples()

# Merge DNAStringSet objects, retaining the first occurrence of each ASV name, and write out
filtered_fastas %>%
  stringr::str_split_1(pattern = " ") %>%
  lapply(Biostrings::readDNAStringSet) %>%
  merge_refseq() %>%
  write_fasta("asvs_filtered.fasta")

if (!identical(phyloseq::rank_names(ps_f), tax_cols)) {
    warning(
        "Taxonomy ranks do not meet staging database expectations.\n",
        "Expected: ", paste(tax_cols, collapse = ", "), "\n",
        "Observed: ", paste(phyloseq::rank_names(ps_f), collapse = ", ")
    )
}

## output merged raw data tables - NOTE: This is memory intensive
data.table::fwrite(
    melt_phyloseq(ps_f),
    "raw_filtered.csv"
)

# Export species level summary of filtered results
data.table::fwrite(
    summarise_phyloseq(ps_f),
    "summary_filtered.csv"
)

## output phyloseq and component data

# Write out phyloseq
saveRDS(ps_f, paste0("ps_filtered.rds"))

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
phyloseq::otu_table(ps_f) %>%
    as("matrix") %>%
    data.table::as.data.table(keep.rownames = "sample_primers") %>%
    data.table::fwrite(file = "seqtab_filtered.csv")

# save taxtab as long tibble (rows = ASV, cols = tax rankings) with unclassified labels converted to NA.
phyloseq::tax_table(ps_f) %>%
    as("matrix") %>%
    unclassified_to_na() %>%
    data.table::as.data.table(keep.rownames = "seq_name") %>%
    data.table::fwrite(file = "taxtab_filtered.csv")

# save samplesheet
samdf <- data.table::as.data.table(
    as(phyloseq::sample_data(ps_f), "data.frame"),
    keep.rownames = "sample_primers"
)
data.table::fwrite(samdf, "samdf_filtered.csv")

## read tracking output from filtered phyloseq
readsout <- {
    ss <- phyloseq::sample_sums(ps_f)
    data.table::data.table(
        sample_primers = names(ss),
        pairs = unname(ss)
    )
}

readsout <- samdf[
    readsout,
    on = "sample_primers"
]

readsout[, stage := "filter_sample_taxon"]

required_cols <- c(
    "stage",
    "sample_primers",
    "read_group",
    "primers",
    "pairs"
)

# Reorder and drop all unwanted metadata columns.
readsout <- readsout[, ..required_cols]

data.table::fwrite(
    readsout,
    "ps_f_readsout.csv"
)

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})