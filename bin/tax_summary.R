#!/usr/bin/env Rscript


## check and define variables
idtaxa <-       readRDS(tax)
idtaxa_ids <-   readRDS(ids)
joint <-        readRDS(joint_file)
# TODO: deal with case where 'joint' does not exist due to BLAST not being performed
# in these cases, set 'joint' to NULL

# TODO: use explicitly defined ranks from IDTAXA database instead of guess
ranks <- c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")

### run R code
OTU_seq <- rownames(idtaxa) # get OTU sequences as vector
OTU_hash <- lapply(OTU_seq, rlang::hash) %>% unlist() # get unique hash for each OTU sequence

idtaxa_summary <- idtaxa_ids %>%
    purrr::map_dfr(function(x){
        taxa <- paste0(x$taxon,"_", x$confidence) # add assignment confidence to taxon name
        taxa[startsWith(taxa, "unclassified_")] <- NA # set unclassified ranks to NA
        data.frame(t(taxa)) %>% # create data frame setting ranks to column names
            magrittr::set_colnames(ranks[1:ncol(.)])
    }) %>%
    mutate_all(function(y){
        name <- y %>% str_remove("_[0-9].*$") # get name of taxon
        conf <- y %>% str_remove("^.*_") %>% # get confidence of taxon, truncated to 5 digits (5 + '.')
            str_trunc(width=6, side="right", ellipsis = "")
        paste0(name, "_", conf, "%") # create new taxon name with confidence in % form
    }) %>%
    mutate_all(~ na_if(., "NA_NA%")) %>% # convert 'NA_NA%' into true NAs
    mutate(
        OTU_seq = OTU_seq, # add OTU sequence as column
        OTU_hash = OTU_hash # add OTU sequence hash as column
    )

if(!is.null(joint)){
    blast_summary <- joint %>% 
        mutate(
            pcr_primers = pcr_primers,
            target_gene = target_gene,
            idtaxa_db = idtaxa_db,
            ref_fasta = ref_fasta,
        ) %>% 
        dplyr::select(
            pcr_primers, 
            target_gene, 
            idtaxa_db, 
            ref_fasta, 
            OTU_seq = OTU, 
            acc,  
            blast_top_hit = blastspp, 
            blast_identity = pident,
            blast_evalue = evalue, 
            blast_total_score = total_score, 
            blast_max_score = max_score, 
            blast_qcov = qcovs
        )
    
    summary_table <- idtaxa_summary %>%
        left_join(blast_summary) %>%
        dplyr::select(any_of(c(
            "OTU_hash",
            "OTU_seq", 
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
    summary_table <- idtaxa_summary %>%
        dplyr::select(any_of(c(
            "OTU_hash",
            "OTU_seq", 
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

out <- paste0(fcid,"_",pcr_primers,"_taxonomic_assignment_summary.csv")
write_csv(summary_table, out)

saveRDS(summary_table, paste0(fcid,"_",pcr_primers,"_taxonomic_assignment_summary.rds"))

# stop(" *** stopped manually *** ") ##########################################
