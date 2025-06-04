#!/usr/bin/env Rscript
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

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "pcr_primers",
    "taxtab_file",
    "seqtab_file",
    "filters_file",
    "fasta_file",
    "loci_params",
    "samdf"
)
lapply(nf_vars, nf_var_check)

## check and define variables
taxtab <- readr::read_csv(taxtab_file)

seqtab_combined <- readr::read_csv(seqtab_file)

filters <- readr::read_csv(filters_file)

seqs_combined <- Biostrings::readDNAStringSet(fasta_file)

samdf <- readr::read_csv(samdf, show_col_types = FALSE)

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

## mutate samdf to add pcr_primers to sample_id, to make consistent with new seqtab format
samdf_renamed <- 
    samdf %>% 
    dplyr::mutate(
        sample_id_orig = sample_id, # save original sample_id as a new column
        sample_id = paste0(sample_id,"_",pcr_primers) # mutate sample_id to add primer id at the end
    )

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
    samdf_renamed %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$sample_id) %>%
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

### outputs -----------------------------------------------------------------------------

# save sequences as .fasta file (with taxonomy in header, in format "seq_name|pcr_primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
seq_names_new <- 
    taxtab %>%
    # ensure sequence name order is same as the DSS object
    dplyr::arrange(factor(seq_name, levels = names(seqs_combined))) %>%
    # add pcr_primers
    dplyr::mutate(pcr_primers = pcr_primers, .after = seq_name) %>%
    # unite columns into a single header string per sequence
    tidyr::unite(col = "lineage", Root:Species, sep = ";") %>%
    tidyr::unite(col = "id", c(seq_name, pcr_primers), sep = "|") %>%
    tidyr::unite(col = "header", c(id, lineage), sep = ";") %>%
    dplyr::pull(header)

seqs_output <- seqs_combined

names(seqs_output) <- seq_names_new

write_fasta(seqs_output, paste0("asvs_unfiltered_", pcr_primers, ".fasta"))  

# save seqtab (filters) as wide tibble (rows = seq_name, columns = sample_id)
readr::write_csv(filters, paste0("filters_",pcr_primers,".csv"))

# save seqtab (data) as wide tibble (rows = seq_name, columns = sample_id)
phyloseq::otu_table(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    readr::write_csv(., paste0("seqtab_unfiltered_",pcr_primers,".csv"))

# save taxtab as tibble (rows = seq_name, columns = taxonomic ranks (Root -> Species))
phyloseq::tax_table(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "seq_name") %>%
    # convert propagated taxonomy to NA values
    dplyr::mutate( 
        dplyr::across(Root:Species, ~ dplyr::if_else(stringr::str_detect(.x, "^\\w__"), NA, .x))
    ) %>%
    readr::write_csv(., paste0("taxtab_unfiltered_",pcr_primers,".csv"))

# save samdf as tibble
phyloseq::sample_data(ps_uf) %>%
    as("matrix") %>%
    tibble::as_tibble() %>%
    readr::write_csv(., paste0("samdf_unfiltered_",pcr_primers,".csv"))

# export raw .csv from phyloseq object
melt_phyloseq(ps_uf) %>%
    readr::write_csv(., paste0("raw_unfiltered_",pcr_primers,".csv"))

# export summary .csv from phyloseq object
summarise_phyloseq(ps_uf) %>%
    dplyr::left_join(., filters, by = c("seq_name", "sequence")) %>%
    dplyr::relocate(seq_name, sequence, Root:Species, chimera_filter, length_filter, phmm_filter, frame_filter) %>%
    readr::write_csv(., paste0("summary_unfiltered_",pcr_primers,".csv"))

# export phyloseq object
saveRDS(ps_uf, paste0("ps_unfiltered_",pcr_primers,".rds"))