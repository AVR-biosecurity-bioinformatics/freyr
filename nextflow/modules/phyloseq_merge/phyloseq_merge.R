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

## check and define variables
ps_unfiltered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_unfiltered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_unfiltered <- lapply(ps_unfiltered, readRDS) # read in phyloseq objects and store as list

# import fasta files
unfiltered_seqs_list <- 
    unfiltered_fastas %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

### run R code

### unfiltered ---------------------------------------------------------------------------------------------

# Combine DNAStringSet objects, retaining the first occurrence of each ASV name.
seqs_uf <- unlist(
    lapply(unfiltered_seqs_list, as.character),
    use.names = TRUE
)
seqs_uf <- Biostrings::DNAStringSet(
    seqs_uf[!duplicated(names(seqs_uf))]
)
write_fasta(seqs_uf, "asvs_unfiltered.fasta")


## merge unfiltered phyloseq objects - These must have disjoint samples
ps_u <- merge_phyloseq_disjoint_samples(ps_unfiltered)

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

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
phyloseq::otu_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "sample_primers") %>%
    readr::write_csv("seqtab_unfiltered.csv")

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
phyloseq::tax_table(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    ) %>%
    readr::write_csv("taxtab_unfiltered.csv")

# save samplesheet
phyloseq::sample_data(ps_u) %>%
    as("matrix") %>%
    tibble::as_tibble() %>%
    readr::write_csv("samdf_unfiltered.csv")

# Clean up
rm(ps_unfiltered)
rm(seqs_uf)
gc()

### filtered ----------------------------------------------------------------------------------------------------------

ps_filtered <- # convert Groovy to R list format
    stringr::str_extract_all(ps_filtered, pattern = "[^\\s,\\[\\]]+") %>% unlist()
ps_filtered <- lapply(ps_filtered, readRDS) # read in phyloseq objects and store as list

filtered_seqs_list <- 
    filtered_fastas %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)


# Combine DNAStringSet objects, retaining the first occurrence of each ASV name.
seqs_f <- unlist(
    lapply(unfiltered_seqs_list, as.character),
    use.names = TRUE
)
seqs_f <- Biostrings::DNAStringSet(
    seqs_f[!duplicated(names(seqs_f))]
)
write_fasta(seqs_f, "asvs_filtered.fasta")

## merge filtered phyloseq objects - These must have disjoint samples
ps_f <- merge_phyloseq_disjoint_samples(ps_filtered)

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
data.table::fwrite(
    data.table::as.data.table(
        as(phyloseq::otu_table(ps_f), "matrix"),
        keep.rownames = "sample_primers"
    ),
    "seqtab_filtered.csv"
)

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
taxtab_out_f <- data.table::as.data.table(
    as(phyloseq::tax_table(ps_f), "matrix"),
    keep.rownames = "seq_name"
)
for (col in tax_cols) {
    taxtab_out_f[
        grepl("^\\w__", get(col)),
        (col) := NA_character_
    ]
}
data.table::fwrite(
    taxtab_out_f,
    "taxtab_filtered.csv"
)

# Check taxonomy table outputs
### TODO: use 'ranks' pipeline parameter (from loci_params?) to set this explicitly rather than guessing
expected_cols <- c("seq_name", tax_cols)
if (!identical(names(taxtab_out_f), expected_cols)) {
    message(
        "Warning: Taxonomy table columns do not meet expectations for the staging database.\n",
        "Database requires the columns: ", paste(expected_cols, collapse = ", ")
    )
}

# save samplesheet
samdf <- data.table::as.data.table(
    as(phyloseq::sample_data(ps_f), "data.frame")
)
data.table::fwrite(samdf, "samdf_filtered.csv")

## read tracking output from filtered phyloseq
sample_sums_f <- phyloseq::sample_sums(ps_f)

readsout <- data.table::data.table(
    sample_primers = names(sample_sums_f),
    pairs = unname(sample_sums_f)
)

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