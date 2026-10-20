#!/usr/bin/env Rscript
.libPaths(.Library)

tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
read_group                  <- args$read_group
fasta                       <- args$fasta
blast_tsv                   <- args$blast_tsv
blast_min_identity          <- args$blast_min_identity
blast_min_coverage          <- args$blast_min_coverage
run_blast                   <- args$run_blast

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "readr",
    "stringr",
    "tibble",
    "tidyr",
    "Biostrings",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

### TODO: Add explicit taxonomic ranks option to loci parameters

### TODO: Do parameter parsing for ranks at an earlier stage, such as the parameter setup module

## check and define variables 

# set variables
run_blast <-            as.logical(run_blast)
blast_min_identity <-   as.numeric(blast_min_identity)
blast_min_coverage <-   as.numeric(blast_min_coverage)

# get sequences
seqmap <- 
    Biostrings::readDNAStringSet(fasta) %>% 
    as.character() %>%
    tibble::enframe(., name = "seq_name", value = "sequence")

# get blast output
blast_out <- 
    readr::read_tsv(
        blast_tsv, 
        col_names = c("qseqid","sseqid","stitle","pident","length","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore","qcovs")
    )

### run code

# check there are some BLAST hits and BLAST is requested
if (nrow(blast_out) > 0 & isTRUE(run_blast)){

    if ( nrow(seqmap) > 0 ) { # if there are ASV sequences, run BLAST

        # get ranks from blast output
        db_rank_count <- blast_out$sseqid %>% stringr::str_count(., ";") %>% unique()

        if (length(db_rank_count) != 1){
        stop(paste0("*** Reference database records appear to have a variable number of taxonomic ranks for primers '",primers,"' ***"))
        }

        if (db_rank_count == 7){
            ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
            message("** BLAST database contains 7 ranks--setting to 'Kingdom>>Species' **")
        } else if (db_rank_count == 8){
            ranks <- c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
            message("** BLAST database contains 8 ranks--setting to 'Root>>Species' **")
        } else if (db_rank_count < 7){
            stop ("*** BLAST database contains fewer than 7 ranks -- please set ranks explicitly using 'params.tax_ranks'! ***")
        } else {
            stop ("*** BLAST database contains more than 8 ranks -- please set ranks explicitly using 'params.tax_ranks'! ***")
        }

        ## create low-stringency BLAST output
        blast_spp_low <- 
            blast_out %>% 
            # convert "!?!?" back to spaces
            dplyr::mutate(
                sseqid = stringr::str_replace_all(sseqid, "\\!\\?\\!\\?", " "),
                stitle = stringr::str_replace_all(stitle, "\\!\\?\\!\\?", " ")
            ) %>% 
            dplyr::filter(!is.na(sseqid)) %>% # Filter out results with no hit
            dplyr::mutate(
                # Normalise coordinates, accounting for reverse-strand alignments.
                q_left  = pmin(qstart, qend),
                q_right = pmax(qstart, qend),
                s_left  = pmin(sstart, send),
                s_right = pmax(sstart, send),
                
                # Count sequence positions spanned by each hit.
                q_align = q_right - q_left + 1,
                s_align = s_right - s_left + 1
            ) %>%
            dplyr::group_by(qseqid, sseqid, stitle) %>%
            dplyr::group_modify(~ {
                # This code loops through hits from best match to worst
                # Top hit is always selected
                # Seperate lower ranked hit(s) are retained only if they dont overlap both query and subject
                # This is common when paired end reads that do not overlap are joined together (i.e. with concat_unmerged in freyr)
                # NOTE: This selection is conservative - a single overlapping base will cause the lower hit to be discarded

                # Sort by bitscore > pident > qcovs
                # Ties are broken reproducibly using coordinates.
                hits <- dplyr::arrange(
                    .x,
                    dplyr::desc(bitscore),
                    dplyr::desc(pident),
                    dplyr::desc(qcovs),
                    q_left, q_right, s_left, s_right
                )
                            
                keep <- integer(nrow(hits))
                n_keep <- 0L
                
                # Loop through hits - Starting with best hit
                for (i in seq_len(nrow(hits))) {
                    # NOTE: top hit is always retianed
                    retained <- keep[seq_len(n_keep)]
                    
                    # Check whether current hit shares any **query** positions with a retained better hit
                    # If it does, it will be removed
                    overlaps_query <- any(
                        hits$q_left[i] <= hits$q_right[retained] &
                        hits$q_right[i] >= hits$q_left[retained]
                    )
                    
                    # Check whether current hit shares any **subject** positions with a retained better hit
                    # If it does, it will be removed
                    overlaps_subject <- any(
                        hits$s_left[i] <= hits$s_right[retained] &
                        hits$s_right[i] >= hits$s_left[retained]
                    )
                    
                    # Retain only hits that overlap on neither sequence.
                    if (!overlaps_query && !overlaps_subject) {
                        n_keep <- n_keep + 1L
                        keep[n_keep] <- i
                    }
                }
                
                hits[keep[seq_len(n_keep)], , drop = FALSE]
            }) %>%
            dplyr::summarise(
                # Adjust each full sequence length for alignment columns containing gaps
                # Use the shorter of the two lengths
                # This handles cases where query is longer than subject. See issue #26
                pident = sum(pident * length) /
                    min(
                        dplyr::first(qlen) + sum(length - q_align),
                        dplyr::first(slen) + sum(length - s_align)
                    ),
                
                # Coverage of retained, non-overlapping query intervals.
                qcovs = 100 * sum(q_align) / dplyr::first(qlen),
                
                max_score = max(bitscore),
                total_score = sum(bitscore),
                evalue = min(evalue),
                .groups = "drop"
            )%>%
            # Initial low stringency filters
            dplyr::filter(pident > 60, qcovs > 80) %>%
            # Retain the best-ranked subject(s) per query, preserving exact ties.
            dplyr::group_by(qseqid) %>%
            dplyr::top_n(1, total_score) %>%
            dplyr::top_n(1, max_score) %>%
            dplyr::top_n(1, qcovs) %>%
            dplyr::top_n(1, pident) %>%
            # each taxonomic rank in its own column
            tidyr::separate(stitle, c("acc", ranks), ";", remove = TRUE) %>%
            dplyr::ungroup()

        # filter by identity and coverage
        blast_spp <- 
            blast_spp_low %>%
            # filter by identity and coverage thresholds
            dplyr::filter(pident >= blast_min_identity, qcovs >= blast_min_coverage) %>% 
            dplyr::group_by(qseqid) %>%
            # add end of species binomial
            dplyr::mutate(
                spp = Species %>% stringr::str_remove("^.* ") %>% stringr::str_remove("^.*_")
            ) %>%
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
    if(!all(blast_spp$seq_name %in% seqmap$seq_name)){
        stop("Number of ASVs classified does not match the number of input ASVs")
    }
        
} else { 

    # if BLAST not requested or there were no hits, produce blast_spp tibble full of NAs
    blast_spp <- 
        seqmap %>%
        dplyr::mutate(
            sequence = NULL,
            blast_genus = NA_character_, 
            blast_spp = NA_character_
        )   

    # NULL output for low stringency BLAST
    blast_spp_low <- NULL

}

# save low stringecy blast output for assignment plot
saveRDS(blast_spp_low, paste0(read_group,"_",primers,"_blast_spp_low.rds"))

# save tibble
readr::write_csv(blast_spp, paste0(read_group,"_",primers,"_blast.csv"))

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})