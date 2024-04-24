#!/usr/bin/env Rscript


## check and define variables



### run R code

stop(" *** stopped manually *** ") ##########################################


## manipulate idtaxa objects (ids and tax)
# idtaxa_ids <- idtaxa$ids # from 'tax_idtaxa'
# idtaxa_tax <- idtaxa$tax # from 'tax_idtaxa'

### NOTE for Alex: 


#### from tar_target(tax_summary)
idtaxa_summary <- tax_idtaxa %>%
    ungroup() %>%
    mutate(
        summary = purrr::map2(idtaxa_ids, idtaxa, ~{
            .x %>%
                purrr::map_dfr(function(x){
                    taxa <- paste0(x$taxon,"_", x$confidence)
                    taxa[startsWith(taxa, "unclassified_")] <- NA
                    data.frame(t(taxa)) %>%
                        magrittr::set_colnames(c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")[1:ncol(.)])
                }) %>%
                mutate_all(function(y){
                    name <- y %>% str_remove("_[0-9].*$")
                    conf <- y %>% str_remove("^.*_") %>%
                        str_trunc(width=6, side="right", ellipsis = "")
                    paste0(name, "_", conf, "%") 
                }) %>%
        mutate_all(~ na_if(., "NA_NA%"))  %>%
        mutate(OTU = rownames(.y))
    })) %>%
    dplyr::select(pcr_primers,idtaxa_db, summary)%>%
    unnest(summary) %>%
    distinct()

if(!any(sapply(assignment_plot$joint, is.null))){
    blast_summary <- assignment_plot %>%
        dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta, joint)%>% 
        unnest(joint)  %>%
        dplyr::select(
            pcr_primers, 
            target_gene, 
            idtaxa_db, 
            ref_fasta, 
            OTU, 
            acc,  
            blast_top_hit = blastspp, 
            blast_identity = pident,
            blast_evalue = evalue, 
            blast_total_score = total_score, 
            blast_max_score = max_score, 
            blast_qcov = qcovs
        )
} else {
    blast_summary <- assignment_plot %>%
        dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta, joint)%>% 
        unnest(joint)
}

summary_table <- idtaxa_summary %>%
    left_join(blast_summary) %>%
    dplyr::select(any_of(
        c("OTU", 
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
        "blast_evalue"
        )
    ))

out <- paste0("output/logs/taxonomic_assignment_summary.csv")
write_csv(summary_table, out)
return(out)