#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dada2",
    "dplyr",
    "readr",
    "S4Vectors",
    "stringr",
    "taxreturn",
    "tibble",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "fcid",
    "pcr_primers",
    "seqtab",
    "target_gene",
    "ref_fasta",
    "blast_min_identity",
    "blast_min_coverage",
    "run_blast"
)
lapply(nf_vars, nf_var_check)

### TODO: Add explicit taxonomic ranks option to loci parameters

### TODO: Do parameter parsing for ranks at an earlier stage, such as the parameter setup module

## check and define variables 
target_gene <-          parse_nf_var_repeat(target_gene)
ref_fasta <-            parse_nf_var_repeat(ref_fasta)
blast_min_identity <-   parse_nf_var_repeat(blast_min_identity)
blast_min_coverage <-   parse_nf_var_repeat(blast_min_coverage)
run_blast <-            parse_nf_var_repeat(run_blast)

quiet <-                FALSE # switch quiet off for now
multithread <-          FALSE # multithreading switched off for now

# extract number of ranks from ref_fasta 
if ( stringr::str_detect(ref_fasta, "\\.fa\\.gz$|\\.fasta\\.gz$") ) { # if compressed
    n_ranks <- readr::read_lines(gzfile(ref_fasta), n_max = 1) %>% # pull first line of .fa.gz
            stringr::str_extract(pattern = ";.*?$") %>% # extract the taxonomic rank information (including first ';')
            stringr::str_count(pattern = "[^;]+") # count the number of ranks between the ';'
} else if ( stringr::str_detect(ref_fasta, "\\.fa$|\\.fasta$") ) { # if uncompressed
    n_ranks <- readr::read_lines(ref_fasta, n_max = 1) %>% # pull first line of .fa.gz
            stringr::str_extract(pattern = ";.*?$") %>% # extract the taxonomic rank information (including first ';')
            stringr::str_count(pattern = "[^;]+") # count the number of ranks between the ';'
} else { # if extension is not expected
    stop ("*** 'ref_fasta' (BLAST database) file must have '.fa(sta)' or '.fa(sta).gz' extension! ***")
}
# set ranks based on number of ranks (only works for 7 or 8 ranks)
if ( n_ranks == 8 ) {
    ranks <- c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    message("** BLAST database contains 8 ranks--setting to 'Root>>Species' **")
} else if ( n_ranks == 7 ) {
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    message("** BLAST database contains 7 ranks--setting to 'Kingdom>>Species' **")
} else if ( n_ranks < 7 ) {
    stop ("*** BLAST database contains fewer than 7 ranks--please set ranks explicitly using 'params.tax_ranks'! ***")
} else {
    stop ("*** BLAST database contains more than 8 ranks--please set ranks explicitly using 'params.tax_ranks'! ***")
}

write(n_ranks, file = paste0(fcid, "_", pcr_primers,"_n_ranks.txt"))

run_blast <-            as.logical(run_blast)
database <-             normalizePath(ref_fasta)
identity <-             as.numeric(blast_min_identity)
coverage <-             as.numeric(blast_min_coverage)
db_name <-              basename(database) %>% stringr::str_remove("_\\.*$")

### run R code
seqtab <- readRDS(seqtab) # read in seqtab


if (isTRUE(run_blast)) { # run BLAST if requested
    
    seqmap <- tibble::enframe(dada2::getSequences(seqtab), name = NULL, value="OTU") %>%
        dplyr::mutate(name = paste0("SV", seq(length(dada2::getSequences(seqtab)))))
    
    seqs <- taxreturn::char2DNAbin(seqmap$OTU)
    names(seqs) <- seqmap$name


    # empty file to debug output
    # EMPTY_FILE <- c()
    
    if ( length(seqs) > 0 ) { # if there are ASV sequences, run BLAST
        ## make low stringency, ident = 60, coverage = 80, then save
        blast_spp_low <- taxreturn::blast_top_hit(
            query = seqs,
            db = database, 
            identity = 60, 
            coverage = 80, 
            evalue = 1e06,
            max_target_seqs = 5, 
            max_hsp = 5, 
            ranks = ranks, 
            delim = ";",
            resolve_ties="all"
            )

        ## save BLAST output for assignment plot
        saveRDS(blast_spp_low, paste0(fcid,"_",pcr_primers,"_blast_spp_low.rds"))

        # filter by identity and coverage
        ### TODO (Alex): update this (you weren't happy with it)
        blast_spp <- blast_spp_low %>% 
            dplyr::filter(pident >= identity, qcovs >= coverage) %>% 
            dplyr::group_by(qseqid) %>%
            dplyr::mutate(spp = Species %>% stringr::str_remove("^.* ")) %>%
            dplyr::reframe(spp = paste(sort(unique(spp)), collapse = "/"), Genus, pident, qcovs, max_score, total_score, evalue) %>%
            dplyr::mutate(binomial = paste(Genus, spp)) %>%
            dplyr::distinct() %>%
            dplyr::group_by(qseqid) %>% # added to resolve issue of returning NAs for Species (add_tally added up all rows ungrouped)
            dplyr::add_tally() %>%
            dplyr::ungroup() %>% 
            dplyr::mutate(binomial =  dplyr::case_when( #Leave unassigned if conflicted at genus level
                n > 1 ~ as.character(NA),
                n == 1 ~ binomial
                )
            ) %>%
            dplyr::select(OTU = qseqid, Genus, Species = binomial, pident, qcovs, max_score, total_score, evalue) %>% 
            dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
            dplyr::filter(!is.na(blast_spp)) 

        # for debug output
        # saveRDS(blast_spp, paste0(fcid,"_",pcr_primers,"_blast_spp.rds"))

        if( nrow(blast_spp) > 0 ) {
        # Transform into taxtab format
        out <- tibble::enframe(dada2::getSequences(seqtab), name=NULL, value="OTU") %>%
            dplyr::left_join(
                blast_spp %>%
                    dplyr::select(name = OTU, Genus = blast_genus, Species = blast_spp) %>%
                    dplyr::left_join(seqmap, by = "name") %>%
                    dplyr::select(-name), 
                by="OTU"
                ) %>%
            tibble::column_to_rownames("OTU") %>%
            as.matrix()
        } else {
            warning(paste0("No Species assigned with BLAST to ", database, " -- have you used the correct database?"))
        out <- tibble::enframe(dada2::getSequences(seqtab), name=NULL, value="OTU") %>%
                dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
                tibble::column_to_rownames("OTU") %>%
                as.matrix()
        }
    } else {
        warning(paste0("No sequences present in seqtab -- BLAST skipped"))
        out <- tibble::enframe(dada2::getSequences(seqtab), name = NULL, value = "OTU") %>%
        dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
        tibble::column_to_rownames("OTU") %>%
        as.matrix()
    }
    
    # Check that output dimensions match input
    if(!all(rownames(out) %in% colnames(seqtab))){
        stop("Number of ASVs classified does not match the number of input ASVs")
    }
    
    # save output
    saveRDS(out, paste0(fcid,"_",pcr_primers,"_",db_name,"_blast.rds"))
    
} else { # if BLAST not requested, produce null output
    
    out <- tibble::enframe(dada2::getSequences(seqtab), name = NULL, value = "OTU") %>%
                    dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
                    tibble::column_to_rownames("OTU") %>%
                    as.matrix()    
    # save output
    saveRDS(out, paste0(fcid,"_",pcr_primers,"_",db_name,"_blast.rds"))

    # save NULL output for assignment plot
    blast_spp_low <- NULL
    saveRDS(blast_spp_low, paste0(fcid,"_",pcr_primers,"_blast_spp_low.rds"))

}

# stop(" *** stopped manually *** ") ##########################################
