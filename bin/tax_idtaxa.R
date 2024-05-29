#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "DECIPHER",
    "dplyr",
    "magrittr",
    "purrr",
    "stringr",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "fcid",
    "pcr_primers",
    "seqtab",
    "idtaxa_confidence",
    "idtaxa_db"
)
lapply(nf_vars, nf_var_check)

## check and define variables 
idtaxa_confidence <-    parse_nf_var_repeat(idtaxa_confidence)
idtaxa_db <-            parse_nf_var_repeat(idtaxa_db)

return_ids <-           TRUE
quiet <-                FALSE # switch quiet off for now
multithread <-          FALSE # multithreading switched off for now
remove_Ns <-            FALSE
ranks <-                c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
### TODO: Check that ranks of IDTAXA and BLAST assignment need to be the same?
### if so, implement similar code to tax_blast.R to extract ranks from database

database <-             idtaxa_db
threshold <-            as.numeric(idtaxa_confidence)


### run R code

seqtab <- readRDS(seqtab) # read in seqtab
trainingSet <- readRDS(normalizePath(database)) # read in training database for IDTAXA

# get the sequences from the seqtab
seqs <- Biostrings::DNAStringSet(dada2::getSequences(seqtab)) # Create a DNAStringSet from the ASVs

if ( length(seqs) > 0 ) { # Stop if seqs are 0

    # Remove any 10bp N bases that were added by concatenating reads reads  
    if(remove_Ns){
        if(any(seqs %>% purrr::map_lgl(~{str_detect(as.character(.x), "NNNNNNNNNN")}))){
        seqs <- Biostrings::DNAStringSet(
            seqs %>% purrr::map_chr(~{str_replace(as.character(.x), "NNNNNNNNNN", "")})
            )
        }
    }

    # Classify 
    ids <- DECIPHER::IdTaxa(
        seqs, 
        trainingSet, 
        processors=1, 
        threshold = threshold, 
        verbose=!quiet, 
        strand = "top"
        ) 

    # Get the filename of that db that we can use to name the output files
    db_name <- basename(database) %>% 
        stringr::str_remove("\\..*$") %>% # remove extension 
        stringr::str_remove("_idtaxa")

    # Check that more than just root has been assigned
    if( any(sapply(ids, function(x){ length(x$taxon) }) > 2)){
        #Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
        tax <- ids %>%
        purrr::map_dfr(function(x){
            taxa <- paste0(x$taxon,"_", x$confidence)
            taxa[startsWith(taxa, "unclassified_")] <- NA
            data.frame(t(taxa)) %>%
            magrittr::set_colnames(ranks[1:ncol(.)])
        }) %>%
        mutate_all(stringr::str_replace,pattern="(?:.(?!_))+$", replacement="") %>%
        magrittr::set_rownames(dada2::getSequences(seqtab))
        # add empty ranks if none were assigned to lower ranks
        tax <- new_bind(tibble::tibble(!!!ranks, .rows = 0, .name_repair = ~ ranks), tax)
    } else {
        warning(paste0("No sequences assigned with IDTAXA to ", database, " have you used the correct database?"))
        tax <- tibble::enframe(dada2::getSequences(seqtab), name=NULL, value="OTU") 
        tax[ranks[1]] <- ranks[1]
        tax[ranks[2:length(ranks)]] <- NA_character_
        tax <- tax %>%
        magrittr::set_rownames(dada2::getSequences(seqtab)) 
    }
} else {
    warning(paste0("No sequences present in seqtab - IDTAXA skipped"))
    tax <- tibble::enframe(dada2::getSequences(seqtab), name=NULL, value="OTU") 
    tax[ranks[1]] <- ranks[1]
    tax[ranks[2:length(ranks)]] <- NA_character_
    tax <- tax %>%
        magrittr::set_rownames(dada2::getSequences(seqtab)) 
}

# Check that output dimensions match input
if(!all(rownames(tax) %in% colnames(seqtab))){
    stop("Number of ASVs classified does not match the number of input ASVs")
}

# Check that all ranks are present
if(!all(colnames(tax) %in% ranks)){
    stop("Number of ranks does not match")
}

# Write out idtaxa objects
saveRDS(tax, paste0(fcid,"_",pcr_primers,"_",db_name,"_idtaxa_tax.rds"))
saveRDS(ids, paste0(fcid,"_",pcr_primers,"_",db_name,"_idtaxa_ids.rds"))

# stop(" *** stopped manually *** ") ##########################################
