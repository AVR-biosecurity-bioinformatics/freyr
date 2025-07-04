#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
read_group                  <- args$read_group
fasta                       <- args$fasta
tax_list                    <- args$tax_list
ids_list                    <- args$ids_list
joint_file                  <- args$joint_file
locus                       <- args$locus
idtaxa_db                   <- args$idtaxa_db
ref_fasta                   <- args$ref_fasta

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "magrittr",
    "purrr",
    "readr",
    "rlang",
    "stringr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables
# read in fasta as tibble
seqs <- 
    Biostrings::readDNAStringSet(fasta) %>% 
    as.character() %>%
    tibble::enframe(name = "seq_name", value = "sequence")

idtaxa_tax <-   
    tax_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) %>% # read in .rds and store as list of tibbles
    dplyr::bind_rows()

idtaxa_ids <-  # lsit of Taxa objects 
    ids_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readRDS)

joint <-        readRDS(joint_file)

# TODO: use explicitly defined ranks from ref fasta database or parameter instead of guess
ranks <- c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")

### run R code
idtaxa_summary <- 
  idtaxa_ids %>%
    # convert each Taxa object into a data frame
    lapply(., function(taxa_obj){
        # transform "Taxa" data frame, one row at a time
        purrr::imap_dfr(.x = taxa_obj, .f = function(x, idx){
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
        })
    }) %>%
    # combine list into one data frame/tibble
    dplyr::bind_rows() %>% 
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
    dplyr::left_join(., seqs, by = "seq_name") %>%
    dplyr::relocate(seq_name, sequence) %>%
     # add metadata
    dplyr::mutate(
        primers = primers,
        locus = locus,
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
            "primers", 
            "locus", 
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
            "primers", 
            "locus", 
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
readr::write_csv(summary_table, paste0(read_group,"_",primers,"_taxonomic_assignment_summary.csv"))

saveRDS(summary_table, paste0(read_group,"_",primers,"_taxonomic_assignment_summary.rds"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})