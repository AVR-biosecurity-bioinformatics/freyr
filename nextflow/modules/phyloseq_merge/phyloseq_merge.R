#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

dataset    <- args$dataset
seqtab_files <- unlist(stringr::str_extract_all(args$seqtab_files, "[^\\s,\\[\\]]+"))
taxtab_files <- unlist(stringr::str_extract_all(args$taxtab_files, "[^\\s,\\[\\]]+"))
samdf_files <- unlist(stringr::str_extract_all(args$samdf_files, "[^\\s,\\[\\]]+" ))
fasta_files <- unlist(stringr::str_extract_all(args$fasta_files, "[^\\s,\\[\\]]+" ))

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
    "phyloseq",
    "stringr",
    "data.table",
    "magrittr",
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

# Read in and merge seqtabs 
# rows = ASVs, columns = samples, as integer matrices
seqtab <- seqtab_files %>%
    lapply(function(file) {
        seqtab_dt <- data.table::fread(file)
        sample_cols <- setdiff(names(seqtab_dt),c("seq_name", "cluster"))
        seqtab <- as.matrix(seqtab_dt[, ..sample_cols])
        rownames(seqtab) <- seqtab_dt$seq_name
        storage.mode(seqtab) <- "integer"
        seqtab
    }) %>%
    merge_seqtabs()

# Read in and merge taxtabs
# rows = ASVs, columns = taxonomy ranks, as character matrices
taxtab <- taxtab_files %>%
    lapply(function(file) {
        taxtab_dt <- data.table::fread(file)
        taxtab <- as.matrix(taxtab_dt[, ..tax_cols])
        rownames(taxtab) <- taxtab_dt$seq_name
        taxtab
    }) %>%
    merge_taxtabs()

# Read in and merge samdfs
# rows = samples, columns = metadata variables, as data frame
samdf <- samdf_files %>%
    lapply(function(file) {
        samdf <- as.data.frame(data.table::fread(file))
        rownames(samdf) <- samdf$sample_primers
        samdf
    }) %>%
    merge_samdfs()

# Read reference sequences and remove any taxonomy information (if present).
seqs <- fasta_files %>%
    lapply(function(file) {
        seqs <- Biostrings::readDNAStringSet(file)
        names(seqs) <- sub("\\|.*$","",names(seqs))
        seqs
    }) %>%
    merge_refseq()

# Write out merged files
data.table::fwrite(
    data.table::as.data.table(
        seqtab,
        keep.rownames = "seq_name"
    ),
    paste0("seqtab_", dataset, ".csv")
)

data.table::fwrite(
    data.table::as.data.table(
        taxtab,
        keep.rownames = "seq_name"
    ),
    paste0("taxtab_", dataset, ".csv")
)
data.table::fwrite(samdf, file = paste0("samdf_", dataset, ".csv"))
write_fasta(seqs, paste0("asvs_", dataset, ".fasta"))

## output merged data tables 
raw_dt <- melt_phyloseq(
    otutab = seqtab,
    taxtab = taxtab,
    samdf = samdf,
    seqs = seqs
)

data.table::fwrite(raw_dt, paste0("raw_", dataset, ".csv"))

rm(raw_dt)
gc()

# Export species level summary of unfiltered results
summary_dt <- summarise_phyloseq(
    otutab = seqtab,
    taxtab = taxtab,
    seqs = seqs
)
data.table::fwrite(summary_dt, paste0("summary_", dataset, ".csv"))
rm(summary_dt)
gc()

## Read-tracking output --------------------------------------------------------

if (dataset == "unfiltered") {
    # Summarise reads classified at each taxonomic rank.
    readsout <- data.table::rbindlist(
        lapply(
            tax_cols,
            function(rank) {

                keep <- !is.na(taxtab[, rank]) &
                    !grepl("^\\w__", taxtab[, rank])

                data.table::data.table(
                    stage = paste0("classified_", tolower(rank)),
                    sample_primers = colnames(seqtab),
                    pairs = colSums(seqtab[keep, , drop = FALSE])
                )
            }
        )
    )

} else {
    # Summarise total reads
    sample_totals <- colSums(seqtab)
    readsout <- data.table::data.table(
        stage = "filter_sample_taxon",
        sample_primers = names(sample_totals),
        pairs = unname(sample_totals)
    )
}

# Add sample metadata required by READ_TRACKING.
readsout <- data.table::as.data.table(samdf)[
    readsout,
    on = "sample_primers"
][
    ,
    .(stage, sample_primers, read_group, primers, pairs)
]


data.table::fwrite(readsout, paste0("ps_", dataset,"_readsout.csv" ))

## create merged phyloseq object
ps_merged <- 
    phyloseq::phyloseq(
        phyloseq::otu_table(seqtab, taxa_are_rows = TRUE), 
        phyloseq::tax_table(taxtab),
        phyloseq::sample_data(samdf),
        phyloseq::refseq(seqs)
    )
saveRDS(ps_merged, paste0("ps_", dataset, ".rds"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})