#!/usr/bin/env Rscript
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

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "phyloseq",
    "purrr",
    "readr",
    "scales",
    "stringr",
    "tibble",
    "data.table",
    "tidyr",
    "vegan",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables
taxtab <- readr::read_csv(taxtab_file)

seqtab_combined <- readr::read_csv(seqtab_file)

filters <- readr::read_csv(filters_file)

seqs_combined <- Biostrings::readDNAStringSet(fasta_file)

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

if(!setequal(seqtab_combined$seq_name, names(seqs_combined))){
    stop("seqtab_combined and seqs_combined do not contain the exact same sequence names")
}


### run R code

## create complete samplesheet with arbitrary metadata from original --samplesheet .csv
samplesheet_complete <- 
    samplesheet_split %>%
    dplyr::left_join(., sample_metadata, by = dplyr::join_by("sample","read_group"))

# create phyloseq-format otu_table (seqtab), ignoring filters -- ASV names are rows
otutab_ps <-   
    seqtab_combined %>%
    dplyr::select(-sequence) %>%
    tibble::column_to_rownames(var = "seq_name") %>%
    as.matrix() %>%
    phyloseq::otu_table(taxa_are_rows = TRUE)

# create dada2-format taxtab
taxtab_ps <-
    taxtab %>%
    tibble::column_to_rownames(var = "seq_name") %>%
    as.matrix() %>%
    phyloseq::tax_table()

# create sample information phyloseq object
samples_ps <-
    samplesheet_complete %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$sample_primers) %>%
    phyloseq::sample_data()

# create refseq object
refseq_ps <- 
    seqs_combined %>%
    phyloseq::refseq()


## create phyloseq object
ps_uf <- 
    phyloseq::phyloseq(
        otutab_ps, 
        taxtab_ps,
        samples_ps,
        refseq_ps
    )

### cluster ASVs ------------------------------------------------------------------------

if (cluster_threshold %>% is.na){
  
    # make clusters "NA" if threshold not given   
    cluster_tibble <- 
        tibble::as_tibble_col(names(seqs_combined), column_name = "seq_name") %>%
        dplyr::mutate(cluster = NA_integer_)

} else {
    
    # cluster at threshold
    set.seed <- 1; clusters <- 
        DECIPHER::Clusterize(
            seqs_combined, 
            cutoff = 1 - (cluster_threshold / 100),
            invertCenters = TRUE,
            # maxPhase1 = 1e5,
            # maxPhase2 = 1e4,
            # maxPhase3 = 1e4, 
            # singleLinkage = FALSE,
            # rareKmers = 100,
            processors = 1
        )

    cluster_tibble <- 
        clusters %>% 
        tibble::as_tibble(rownames = "seq_name") %>%
        dplyr::mutate(cluster = abs(cluster))

}

### outputs -----------------------------------------------------------------------------

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seq_names_new <- 
    taxtab %>%
    # ensure sequence name order is same as the DSS object
    dplyr::arrange(factor(seq_name, levels = names(seqs_combined))) %>%
    # add primers
    dplyr::mutate(primers = primers, .after = seq_name) %>%
    # unite columns into a single header string per sequence
    tidyr::unite(col = "lineage", Root:Species, sep = ";") %>%
    tidyr::unite(col = "id", c(seq_name, primers), sep = "|") %>%
    tidyr::unite(col = "header", c(id, lineage), sep = ";") %>%
    dplyr::pull(header)

seqs_output <- seqs_combined

names(seqs_output) <- seq_names_new

write_fasta(seqs_output, paste0("asvs_unfiltered_", primers, ".fasta"))  

# save seqtab (filters) as wide tibble (rows = seq_name, columns = sample_primers)
readr::write_csv(filters, paste0("filters_",primers,".csv"))

# save seqtab (data) as wide tibble (rows = seq_name, columns = cluster, all of sample_primers)
phyloseq::otu_table(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = seq_name) %>%
    readr::write_csv(., paste0("seqtab_unfiltered_",primers,".csv"))

# save taxtab as tibble (rows = seq_name, columns = taxonomic ranks (Root -> Species))
phyloseq::tax_table(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = seq_name) %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    ) %>%
    readr::write_csv(., paste0("taxtab_unfiltered_",primers,".csv"))

# save samdf as tibble
phyloseq::sample_data(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble() %>%
    readr::write_csv(., paste0("samdf_unfiltered_",primers,".csv"))

# export raw .csv from phyloseq object
melt_phyloseq(ps_uf) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(cluster, .after = sequence) %>%
    readr::write_csv(., paste0("raw_unfiltered_",primers,".csv"))

# export summary .csv from phyloseq object
summarise_phyloseq(ps_uf) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(., filters, by = c("seq_name", "sequence")) %>%
    dplyr::left_join(., cluster_tibble, by = "seq_name") %>% 
    dplyr::relocate(seq_name, sequence, cluster, Root:Species, chimera_filter, length_filter, phmm_filter, frame_filter) %>%
    readr::write_csv(., paste0("summary_unfiltered_",primers,".csv"))

# export phyloseq object
saveRDS(ps_uf, paste0("ps_unfiltered_",primers,".rds"))

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})