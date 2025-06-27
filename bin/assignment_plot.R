#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "readr",
    "seqateurs",
    "taxreturn",
    "tibble",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "read_group",     
    "primers",   
    "fasta",
    "blast_list",  
    "tax",           
    "idtaxa_db",  
    "ref_fasta"  
)
lapply(nf_vars, nf_var_check)

## read in files
fasta <-   Biostrings::readDNAStringSet(fasta)

blast_input <- 
    blast_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readRDS) %>% # read in .rds and store as list of tibbles
    dplyr::bind_rows()

tax_input <-      readr::read_csv(tax)

### run R code

# if blast output isn't NULL...
if ( nrow(blast_input) > 0 ){

    # convert blast 'resolve_ties="all"' output to 'resolve_ties="first"'
    blast <- blast_input %>%
        dplyr::group_by(qseqid) %>% 
        dplyr::mutate(row_n = dplyr::row_number()) %>%
        dplyr::top_n(1, row_n) %>% # Break ties by position
        dplyr::select(-row_n) %>%
        dplyr::ungroup()

    # add "lowest" (lowest assigned rank) to tax tibble
    tax <- 
      tax_input %>%
        dplyr::mutate(lowest = dplyr::case_when(
          stringr::str_detect(Kingdom, "__")  ~ "Root",
          stringr::str_detect(Phylum, "__")  ~ "Kingdom",
          stringr::str_detect(Class, "__")  ~ "Phylum",
          stringr::str_detect(Order, "__")  ~ "Class",
          stringr::str_detect(Family, "__")  ~ "Order",
          stringr::str_detect(Genus, "__")  ~ "Family",
          stringr::str_detect(Species, "__")  ~ "Genus",
          .default = "Species"
          )
        )
        
    # make seqmap 
    seqmap <- tibble::enframe(fasta %>% as.character(), name = "seq_name", value="sequence") 
    # get ASV sequences as named vector
    seqs <- fasta %>% as.character()

    # add sequences to blast output and rename columns
    blast <- blast %>% 
        dplyr::mutate(blastspp = Species) %>%
        dplyr::select(seq_name = qseqid, acc, blastspp, pident, total_score, max_score, evalue, qcovs) %>%
        dplyr::left_join(., seqmap, by = dplyr::join_by("seq_name")) 

    # combine blast output and tax table
    if ( nrow(blast) > 0 & nrow(tax) > 0 ) {
        joint <- blast %>% 
            dplyr::left_join(tax, by = "seq_name")

    } else {
        joint <- NULL
    }
} else {
    joint <- NULL
}

# save joint object
saveRDS(joint, paste0(read_group,"_",primers,"_joint.rds"))

# make assignment plot
if ( !is.null(joint) ) {
    # colours for plot
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

    plot <- joint %>%
        dplyr::select(pident, rank = lowest) %>%
        dplyr::mutate(rank = factor(rank, levels = c("Root","Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
        ggplot(aes(x=pident, fill=rank)) + 
        geom_histogram(colour="black", binwidth = 1, position = "stack") + 
        labs(
            title = paste0(read_group, "  ", primers, " Top hit identity distribution"),
            subtitle = paste0("IDTAXA database:", idtaxa_db, "\nBLAST database:", ref_fasta),
            x = "BLAST top hit % identity",
            y = "ASVs"
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
            panel.grid = element_line(size = rel(0.5))
        ) 
} else { plot <- NULL }

# write plot
plot_filename <- paste0(read_group,"_",primers,"_taxonomic_assignment_summary.pdf")
if ( !is.null(plot) ){
    pdf(file = plot_filename, width = 11, height = 8 , paper="a4r")
    print(plot)
    try(dev.off(), silent=TRUE)
} else {
    pdf(file = plot_filename, width = 11, height = 8 , paper="a4r")
    plot.new()
    text(x=.5, y=.5, "No blast hits to reference fasta -- assignment plot not created") 
    try(dev.off(), silent=TRUE)
}
