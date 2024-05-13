#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    # "Biostrings",
    # "bs4Dash",
    # "clustermq",
    # "dada2",
    # "DECIPHER",
    "dplyr",
    # "future",
    # "ggplot2",
    # "gridExtra",
    # "gt",
    # "magrittr",
    # "markdown",
    # "ngsReports",
    # "patchwork",
    # "phyloseq",
    # "pingr",
    "purrr",
    # "readr",
    # "rlang",
    # "rstudioapi",
    # "savR",
    # "scales",
    # "seqateurs",
    # "shiny",
    # "shinybusy",
    # "shinyWidgets",
    # "ShortRead",
    "stringr",
    # "taxreturn",
    "tibble",
    # "tidyr",
    # "vegan",
    # "visNetwork",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables 
target_gene <-          parse_nf_var_repeat(target_gene)

propagate_tax <-        TRUE

### run R code      # derived from step_join_tax_blast() in functions.R and tar_target(joint_tax)
idtaxa_output <-              readRDS(idtaxa_output)
blast_output <-               readRDS(blast_output)
seqtab <-                     readRDS(seqtab)

# merge idtaxa outputs
if ( length(idtaxa_output) == 1 ) {
    idtaxa_output <- idtaxa_output[[1]]
} else if( length(idtaxa_output) == 2 ) {
    idtaxa_output <- coalesce_tax(idtaxa_output[[1]], idtaxa_output[[2]])
} else if( length(idtaxa_output) == 3 ) {
    temptax <- coalesce_tax(idtaxa_output[[1]], idtaxa_output[[2]])
    idtaxa_output <- coalesce_tax(temptax, idtaxa_output[[3]])
}

# Check that idtaxa_output dimensions match input
if(!all(rownames(idtaxa_output) %in% colnames(seqtab))){
    stop("Number of ASVs classified does not match the number of input ASVs [idtaxa_output]")
}

# merge blast outputs
if( length(blast_output) == 1 ) {
    blast_output <- blast_output[[1]]
} else if( length(blast_output) == 2 ) {
    blast_output <- coalesce_tax(blast_output[[1]], blast_output[[2]])
} else if( length(blast_output) == 3 ) {
    temptax <- coalesce_tax(blast_output[[1]], blast_output[[2]])
    blast_output <- coalesce_tax(temptax, blast_output[[3]])
}

# Check that blast_output dimensions match input
if(!all(rownames(blast_output) %in% colnames(seqtab))){
    stop("Number of ASVs classified does not match the number of input ASVs [blast_output]")
}

# another check
if(!all(rownames(idtaxa_output) %in% rownames(blast_output))){
    stop("ASVs in idtaxa_output and blast_output do not match")
}

# try to merge IDTAXA and BLAST assignments
if(nrow(idtaxa_output) > 0 & nrow(blast_output) > 0){
    #Set BLAST ids where idtaxa_output isn't at genus level to NA
    disagreements <- !blast_output[,1] == idtaxa_output[,7]
    disagreements[is.na(disagreements)] <- TRUE
    blast_output[,1][disagreements] <- NA
    blast_output[,2][disagreements] <- NA
    
    # Join the two taxonomies, prefering names from the tax
    tax_blast <- coalesce_tax(idtaxa_output, blast_output)

    if ( propagate_tax ) {
        tax_blast <- tax_blast %>%
            seqateurs::na_to_unclassified() # Propagate high order ranks to unassigned ASVs
    }
} else {
    warning("Either idtaxa_output or blast_output is empty")
    tax_blast <- as.matrix(idtaxa_output)
}

# another check
if ( !all(rownames(idtaxa_output) %in% rownames(tax_blast)) ) {
    stop("Number of ASVs output does not match the number of input ASVs")
}

# Write taxonomy table for that db to disk
saveRDS(tax_blast, paste0(fcid,"_",pcr_primers,"_taxblast.rds"))

