#!/usr/bin/env Rscript


## check and define variables


### run R code

## for assignment plot, need:
# - 'filtered_seqtab', from FILTER_SEQTAB
# - 'tax', the taxblast object from JOINT_TAX, with additional modifications
# - 'blast', the blast output object from TAX_BLAST with additional piped modifications
# - 



### add resolve_ties= 'first' code filtering to the blast_spp_low output
# blast_spp_low <- blast_spp_low %>%
#     dplyr::mutate(row_n = dplyr::row_number()) %>%
#     dplyr::top_n(1, row_n) %>% # Break ties by position
#     dplyr::select(-row_n) %>%
#     dplyr::ungroup()


stop(" *** stopped manually *** ") ##########################################


## from _targets 
temp_samdf3_grouped %>%
    dplyr::select(-one_of("target_gene", "idtaxa_db"))%>%
    dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta)) %>%
    tidyr::separate_rows(ref_fasta, sep=";") %>%
    dplyr::group_by(fcid, target_gene, pcr_primers, idtaxa_db, ref_fasta) %>%
    tidyr::nest()  %>% 
    dplyr::mutate(ref_fasta2 = purrr::map(ref_fasta, ~{
        ref_fasta_tracked[stringr::str_detect(ref_fasta_tracked, .x)]
    }))  %>%
    unnest(ref_fasta2) %>%
    dplyr::mutate(
        filtered_seqtab = purrr::map2(
            fcid, pcr_primers, 
            ~{
                #TODO this may need to be merged by PCR primer
                readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path,  paste0(.x,"_", .y, "_seqtab.cleaned.rds"))]) 
    }))%>% 
    dplyr::mutate(
        tax = purrr::map2(
            fcid,pcr_primers, 
            ~{
                readRDS(joint_tax[stringr::str_detect(joint_tax, paste0(.x,"_",.y,"_taxblast.rds"))])%>% 
                seqateurs::unclassified_to_na(rownames=FALSE) %>%
                dplyr::mutate(lowest = seqateurs::lowest_classified(.)) 
    })) %>%
    dplyr::mutate(
        blast = purrr::pmap(
            list(pcr_primers, filtered_seqtab, ref_fasta2),
            .f = ~{
                seqmap <- tibble::enframe(getSequences(..2), name = NULL, value="OTU") %>%
                    mutate(name = paste0("SV", seq(length(getSequences(..2)))))
                seqs <- taxreturn::char2DNAbin(seqmap$OTU)
                names(seqs) <- seqmap$name
                



                if(length(seqs) > 0) {
                    blast_top_hit(
                        query = seqs,
                        db = ..3,
                        identity=60,
                        coverage=80,
                        ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")
                        ) %>% 
                    mutate(blastspp = paste0(Genus, " ", Species)) %>%
                    dplyr::select(name = qseqid, acc, blastspp, pident, total_score, max_score, evalue, qcovs) %>%
                    left_join(seqmap) %>%
                    dplyr::select(-name)
                } else {
                    tibble::enframe(seqs, name=NULL, value="OTU") %>%
                    dplyr::mutate(
                        acc = NA_character_,
                        blastspp = NA_character_,
                        pident = NA_real_, 
                        length = NA_integer_,
                        evalue = NA_real_,
                        qcovs = NA_integer_
                        ) 
                }
    })) %>%
    dplyr::mutate(
        joint = purrr::pmap(
            list(blast, tax),
            .f = ~{
                if(nrow(..1) > 0 & nrow(..2) > 0){
                ..1 %>%
                    dplyr::left_join(..2, by="OTU")
                } else {
                    NULL
                }
    })) %>%
    dplyr::mutate(
        plot = purrr::pmap(
            list(fcid, pcr_primers, joint, idtaxa_db, ref_fasta),
            .f= ~{
                # ADD TITLES HERE!
                if(!is.null(..3)){
                    cols <- c(
                        Root = "#D53E4F",
                        Kingdom = "#F46D43",
                        Phylum = "#FDAE61",
                        Class = "#FEE08B",
                        Order = "#E6F598",
                        Family = "#ABDDA4",
                        Genus = "#66C2A5",
                        Species = "#3288BD"
                        ) 
                    ..3 %>%
                        dplyr::select(pident, rank = lowest) %>%
                        dplyr::mutate(rank = factor(rank, levels = c("Root","Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
                        ggplot(aes(x=pident, fill=rank))+ 
                        geom_histogram(colour="black", binwidth = 1, position = "stack") + 
                        labs(
                            title = paste0(..1, "  ", ..2, " Top hit identity distribution"),
                            subtitle = paste0("IDTAXA database:", ..4, " BLAST database:", ..5),
                            x = "BLAST top hit % identity",
                            y = "Sequence Variants"
                            ) + 
                        scale_x_continuous(breaks=seq(60,100,2)) +
                        scale_fill_manual(name = "Taxonomic \nAssignment", values = cols)+
                        theme_bw()+
                        theme(
                        strip.background = element_rect(colour = "black", fill = "lightgray"),
                        strip.text = element_text(size=9, family = ""),
                        plot.background = element_blank(),
                        text = element_text(size=9, family = ""),
                        axis.text = element_text(size=8, family = ""),
                        legend.position = "right",
                        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                        panel.grid = element_line(size = rel(0.5)),
                        ) 
                    } else {
                        NULL
                    }
                }))
