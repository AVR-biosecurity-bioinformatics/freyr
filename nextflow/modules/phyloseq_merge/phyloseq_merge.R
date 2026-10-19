#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

dataset    <- args$dataset
ps_files   <- args$ps_files
fasta_files <- args$fasta_files

if (!dataset %in% c("filtered", "unfiltered")) {
    stop("dataset must be 'filtered' or 'unfiltered'; received: ",dataset)
}

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

# Define tax_cols
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

## merge unfiltered phyloseq objects - These must have disjoint samples
ps_merged <- ps_files %>%
  stringr::str_extract_all(pattern = "[^\\s,\\[\\]]+") %>%
  unlist() %>%
  lapply(readRDS) %>%
  merge_phyloseq_disjoint_samples()

# Merge DNAStringSet objects, retaining the first occurrence of each ASV name, and write out
fasta_files %>%
  stringr::str_split_1(pattern = " ") %>%
  lapply(Biostrings::readDNAStringSet) %>%
  merge_refseq() %>%
  write_fasta(paste0("asvs_", dataset, ".fasta"))

## output merged data tables 
data.table::fwrite(
    melt_phyloseq(ps_merged),
    paste0("raw_", dataset, ".csv")
)

# Export species level summary of unfiltered results
ps_merged_summary <- summarise_phyloseq(ps_merged)
data.table::fwrite(ps_merged_summary, paste0("summary_", dataset, ".csv"))

## output phyloseq and component data

# Write out whole phyloseq object
saveRDS(ps_merged, paste0("ps_", dataset, ".rds"))

## Read-tracking output --------------------------------------------------------

# Convert sample metadata without duplicating sample_primers.
samdf_df <- as(phyloseq::sample_data(ps_merged), "data.frame")
samdf <- data.table::as.data.table(samdf_df)

if (dataset == "unfiltered") {

    # Summarise reads classified at each taxonomic rank.
    sample_cols <- samdf$sample_primers

    abundance <- as.matrix(ps_merged_summary[ , ..sample_cols])
    storage.mode(abundance) <- "integer"

    readsout <- data.table::rbindlist(
        lapply(
            tax_cols,
            function(rank) {

                classified <- (
                    !is.na(ps_merged_summary[[rank]]) &
                    !grepl(
                        "^\\w__",
                        ps_merged_summary[[rank]]
                    )
                )

                data.table::data.table(
                    stage = paste0(
                        "classified_",
                        tolower(rank)
                    ),
                    sample_primers = sample_cols,
                    pairs = colSums(
                        abundance[
                            classified,
                            ,
                            drop = FALSE
                        ]
                    )
                )
            }
        )
    )

} else {

    # Summarise reads remaining after all taxon-, observation- and
    # sample-level filtering.
    sample_totals <- phyloseq::sample_sums(ps_merged)

    readsout <- data.table::data.table(
        stage = "filter_sample_taxon",
        sample_primers = names(sample_totals),
        pairs = unname(sample_totals)
    )
}

# Add sample metadata required by READ_TRACKING.
readsout <- samdf[readsout, on = "sample_primers"][
    ,
    .(
        stage,
        sample_primers,
        read_group,
        primers,
        pairs
    )
]

data.table::fwrite(
    readsout,
    paste0(
        "ps_",
        dataset,
        "_readsout.csv"
    )
)

# Clean up large summary objects.
rm(readsout, ps_merged_summary)
if (exists("abundance")) {
    # Only exists in unfiltered ouput
    rm(abundance)
}
gc()

# save seqtab as wide tibble (rows = sample_primers, cols = ASV name (hash), cells = abundance)
phyloseq::otu_table(ps_merged) %>%
    as("matrix") %>%
    data.table::as.data.table(
        keep.rownames = "sample_primers"
    ) %>%
    data.table::fwrite(file = paste0("seqtab_", dataset, ".csv"))

# save taxtab as long tibble (rows = ASV, cols = tax rankings)
phyloseq::tax_table(ps_merged) %>%
    as("matrix") %>%
    unclassified_to_na() %>%
    data.table::as.data.table(keep.rownames = "seq_name") %>%
    data.table::fwrite(file = paste0("taxtab_", dataset, ".csv"))

# save samplesheet
phyloseq::sample_data(ps_merged) %>%
    as("data.frame") %>%
    data.table::as.data.table(keep.rownames = "sample_primers") %>%
    data.table::fwrite(file = paste0("samdf_", dataset, ".csv"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})