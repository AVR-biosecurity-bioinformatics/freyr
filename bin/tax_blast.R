#!/usr/bin/env Rscript

## check and define variables 
target_gene <-          parse_nf_var_repeat(target_gene)
ref_fasta <-            parse_nf_var_repeat(ref_fasta)
blast_min_identity <-   parse_nf_var_repeat(blast_min_identity)
blast_min_coverage <-   parse_nf_var_repeat(blast_min_coverage)
run_blast <-            parse_nf_var_repeat(run_blast)

quiet <-                FALSE # switch quiet off for now
multithread <-          FALSE # multithreading switched off for now
ranks <-                c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

run_blast <-            as.logical(run_blast)
database <-             ref_fasta
identity <-             blast_min_identity
coverage <-             blast_min_coverage

print(ref_fasta)
print(target_gene)

stop(" *** stopped manually *** ") ##########################################
### run R code
seqtab <- readRDS(seqtab) # read in seqtab


if (isTRUE(run_blast)) { # run BLAST if requested
    
    seqmap <- tibble::enframe(getSequences(seqtab), name = NULL, value="OTU") %>%
        mutate(name = paste0("SV", seq(length(getSequences(seqtab)))))
    
    seqs <- taxreturn::char2DNAbin(seqmap$OTU)
    names(seqs) <- seqmap$name
    
    # Get the filename of that db that we can use to name the output files
    db_name <- basename(database) %>% stringr::str_remove("_.*$")
    
    if ( length(seqs) > 0 ) { # if there are ASV sequences, run BLAST
        
        blast_spp <- taxreturn::blast_assign_species(
            query = seqs,
            db = database, 
            identity = identity, 
            coverage = coverage, 
            evalue = 1e06,
            max_target_seqs = 5, 
            max_hsp = 5, 
            ranks = ranks, 
            delim = ";"
            ) %>%
            dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
            dplyr::filter(!is.na(blast_spp)) 
        
        if( nrow(blast_spp) > 0 ) {
        # Transform into taxtab format
        out <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") %>%
            dplyr::left_join(
                blast_spp %>%
                dplyr::select(name = OTU, Genus = blast_genus, Species = blast_spp) %>%
                left_join(seqmap) %>%
                dplyr::select(-name), by="OTU"
                ) %>%
            column_to_rownames("OTU") %>%
            as.matrix()
        } else {
            warning(paste0("No Species assigned with BLAST to ", database, " -- have you used the correct database?"))
        out <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") %>%
                dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
                column_to_rownames("OTU") %>%
                as.matrix()
        }
    } else {
        warning(paste0("No sequences present in seqtab -- BLAST skipped"))
        out <- tibble::enframe(getSequences(seqtab), name = NULL, value = "OTU") %>%
        dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
        column_to_rownames("OTU") %>%
        as.matrix()
    }
    
    # Check that output dimensions match input
    if(!all(rownames(out) %in% colnames(seqtab))){
        stop("Number of ASVs classified does not match the number of input ASVs")
    }
    
    if(!is.null(output)){
        saveRDS(out, paste0(fcid,"_",pcr_primers,"_",db_name,"_blast.rds"))
        }
    
} else { # if BLAST not requested, produce null output
    
    out <- tibble::enframe(getSequences(seqtab), name = NULL, value = "OTU") %>%
                    dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
                    column_to_rownames("OTU") %>%
                    as.matrix()
    
    saveRDS(out, paste0(fcid,"_",pcr_primers,"_",db_name,"_blast.rds"))

}

# stop(" *** stopped manually *** ") ##########################################
