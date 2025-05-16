#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "purrr",
    "readr",
    "rlang",
    "stringr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "pcr_primers",
    "fcid",
    "tax_file",
    "ids_file",
    "joint_file",
    "target_gene",
    "idtaxa_db",
    "ref_fasta"
)
lapply(nf_vars, nf_var_check)

## check and define variables
idtaxa_tax <-   readr::read_csv(tax_file)
idtaxa_ids <-   readRDS(ids_file)
joint <-        readRDS(joint_file)

# TODO: use explicitly defined ranks from ref fasta database or parameter instead of guess
ranks <- c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")

### run R code
idtaxa_summary <- 
  idtaxa_ids %>%
    # transform "Taxa" data frame, one row at a time
    purrr::imap_dfr(function(x, idx){
        # get assigned taxa as vector
        taxa <- paste0(x$taxon, "_", x$confidence)
        # make unclassified taxa NA
        taxa[startsWith(taxa, "unclassified_")] <- NA
        # make a data frame with the ranks as columns and taxa as values, seq_name as column too
        out_df <- 
            data.frame(t(taxa)) %>% 
            magrittr::set_colnames(ranks[1:ncol(.)]) %>%
            dplyr::mutate(seq_name = idx) %>%
            dplyr::relocate(seq_name)

        return(out_df)
    }) %>%
    dplyr::mutate(
      # add confidence as % to end of taxon name
      dplyr::across(Root:Species, ~{
        tax_name <- .x %>% stringr::str_remove("_[0-9]+.*$")
        tax_conf <- .x %>% stringr::str_remove(paste0(tax_name,"_")) %>% as.numeric() %>% round(., digits = 1)
        return(paste0(tax_name, "__", tax_conf, "%"))
      }),
      # convert "NA__NA%" into true NA
      dplyr::across(Root:Species, ~ dplyr::na_if(.x, "NA__NA%"))
    ) %>%
    # add ASV sequence
    dplyr::left_join(., joint %>% dplyr::select(seq_name, sequence), by = "seq_name") %>%
    dplyr::relocate(seq_name, sequence) %>%
     # add metadata
    dplyr::mutate(
        pcr_primers = pcr_primers,
        target_gene = target_gene,
        idtaxa_db = idtaxa_db,
        ref_fasta = ref_fasta,
    )

if(!is.null(joint)){
    blast_summary <- 
      joint %>% 
        dplyr::select(
            seq_name,
            acc,  
            blast_top_hit = blastspp, 
            blast_identity = pident,
            blast_evalue = evalue, 
            blast_total_score = total_score, 
            blast_max_score = max_score, 
            blast_qcov = qcovs
        )

    summary_table <- 
      idtaxa_summary %>%
        dplyr::left_join(., blast_summary, by = "seq_name") %>%
        dplyr::select(tidyselect::any_of(c(
            "seq_name",
            "sequence", 
            "pcr_primers", 
            "target_gene", 
            "idtaxa_db", 
            "ref_fasta", 
            "Root",
            "Kingdom", 
            "Phylum",
            "Class", 
            "Order", 
            "Family",
            "Genus",
            "Species",
            "blast_top_hit",
            "blast_identity",
            "blast_qcov",
            "blast_evalue",
            "acc"
            )
        ))

} else {
    summary_table <- 
      idtaxa_summary %>%
        dplyr::select(tidyselect::any_of(c(
            "seq_name",
            "sequence", 
            "pcr_primers", 
            "target_gene", 
            "idtaxa_db", 
            "ref_fasta", 
            "Root",
            "Kingdom", 
            "Phylum",
            "Class", 
            "Order", 
            "Family",
            "Genus",
            "Species",
            "blast_top_hit",
            "blast_identity",
            "blast_qcov",
            "blast_evalue",
            "acc"
            )
        ))
}

# write out
readr::write_csv(summary_table, paste0(fcid,"_",pcr_primers,"_taxonomic_assignment_summary.csv"))

saveRDS(summary_table, paste0(fcid,"_",pcr_primers,"_taxonomic_assignment_summary.rds"))

# stop(" *** stopped manually *** ") ##########################################
