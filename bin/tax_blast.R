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
    "fasta",
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
seqs <-  Biostrings::readDNAStringSet(fasta) %>% as.character() # read in fasta


if (isTRUE(run_blast)) { # run BLAST if requested
    
    seqmap <- seqs %>% tibble::enframe(., name = "seq_name", value = "sequence")
    
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
        blast_spp <- 
            blast_spp_low %>%
            # filter by identity and coverage thresholds
            dplyr::filter(pident >= identity, qcovs >= coverage) %>% 
            dplyr::group_by(qseqid) %>%
            # add end of species binomial
            dplyr::mutate(spp = Species %>%
                            stringr::str_remove("^.* ") %>%
                            stringr::str_remove("^.*_")) %>%
            # create "/" for species name when multiple are best hits for one species, remove old Species name
            dplyr::reframe(spp = paste(sort(unique(spp)), collapse = "/"), Genus, pident, qcovs, max_score, total_score, evalue) %>%
            # create new binomial if its not present already
            dplyr::mutate(binomial = paste(Genus, spp)) %>%
            # remove duplicate IDs for the same sequence
            dplyr::distinct() %>%
            dplyr::group_by(qseqid) %>% # added to resolve issue of returning NAs for Species (add_tally added up all rows ungrouped)
            # count number of best hits
            dplyr::add_tally() %>%
            dplyr::ungroup() %>% 
            # make sure binomial is NA if more than one Genus is assigned to a species
            dplyr::mutate(
                binomial =  dplyr::case_when( #Leave unassigned if conflicted at genus level
                    n > 1 ~ as.character(NA),
                    n == 1 ~ binomial
                )
            ) %>%
            # remove unwanted columns, making new Species column modified binomial
            dplyr::select(seq_name = qseqid, Genus, Species = binomial, pident, qcovs, max_score, total_score, evalue) %>% 
            # rename assignment columns
            dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
            # remove sequences without assignment to species level
            dplyr::filter(!is.na(blast_spp)) 

    } else {
        stop("No sequences present in input FASTA")
    }
    
    # Check that output sequences match input
    if(!all(blast_spp$seq_name %in% names(seqs))){
        stop("Number of ASVs classified does not match the number of input ASVs")
    }
        
} else { 

    # if BLAST not requested, produce blast_spp tibble full of NAs
    blast_spp <- 
        seqmap %>%
        dplyr::mutate(
            sequence = NULL,
            blast_genus = NA_character_, 
            blast_spp = NA_character_
        )   

    # save NULL output for assignment plot
    blast_spp_low <- NULL
    saveRDS(blast_spp_low, paste0(fcid,"_",pcr_primers,"_blast_spp_low.rds"))

}

# save tibble
readr::write_csv(blast_spp, paste0(fcid,"_",pcr_primers,"_",db_name,"_blast.csv"))

# stop(" *** stopped manually *** ") ##########################################
