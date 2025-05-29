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
    "seqtab_list",
    "fasta_list",
    "loci_params",
    "samdf"
)
lapply(nf_vars, nf_var_check)

## check and define variables
taxtab <- readr::read_csv(taxtab_file)

seqtab_list <- 
    seqtab_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) # read in seqtabs and store as list of tibbles

fasta_list <- 
    fasta_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

samdf <- readr::read_csv(samdf, show_col_types = FALSE)

### run R code


# combine seqtab tibbles
seqtab_combined <- 
    seqtab_list %>%
    # pivot each tibble longer
    lapply(
        .,
        function(x){ # per tibble
            x %>%
                tidyr::pivot_longer(
                cols = !c(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter),
                names_to = "sample_id",
                values_to = "abundance"
                )
        }
    ) %>%
    # bind tibbles together now columns all match
    dplyr::bind_rows() %>%
    # pivot wider, filling missing abundance with 0
    tidyr::pivot_wider(
        names_from = sample_id,
        values_from = abundance, 
        values_fill = 0
    )

# combine fasta DSS objects, removing redundant sequences
seqs_combined <- 
    fasta_list %>%
    lapply(., as.character) %>%
    unlist(use.names = T) %>% 
    tibble::enframe() %>%
    dplyr::group_by(name, value) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet()

## check taxonomy, seqtab and fasta all have the same sequences, none are missing
if(!setequal(taxtab$seq_name, seqtab_combined$seq_name)){
    stop("taxtab and seqtab_combined do not contain the exact same sequence names")
}
if(!setequal(taxtab$seq_name, names(seqs_combined))){
    stop("taxtab and seqs_combined do not contain the exact same sequence names")
}
if(!setequal(names(seqs_combined), seqtab_combined$seq_name)){
    stop("seqs_combined and seqtab_combined do not contain the exact same sequence names")
}

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
    dplyr::select(-c(sequence, chimera_filter, length_filter, phmm_filter, frame_filter)) %>%
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
seqtab_combined %>%
    dplyr::select(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter) %>%
    readr::write_csv(., paste0("filters_",pcr_primers,".csv"))

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
    dplyr::left_join(
        .,
        seqtab_combined %>% 
            dplyr::select(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter),
        by = c("seq_name", "sequence")
    ) %>%
    dplyr::relocate(seq_name, sequence, Root:Species, chimera_filter:frame_filter) %>%
    readr::write_csv(., paste0("summary_unfiltered_",pcr_primers,".csv"))

# export phyloseq object
saveRDS(ps_uf, paste0("ps_unfiltered_",pcr_primers,".rds"))

# stop(" *** stopped manually *** ") ##########################################
