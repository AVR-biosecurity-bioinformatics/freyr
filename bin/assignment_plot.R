#!/usr/bin/env Rscript


## check and define variables
seqtab <-   readRDS(seqtab)
blast <-    readRDS(blast) # blast is the low-stringency blast output from TAX_BLAST
tax <-      readRDS(tax)

# convert blast 'resolve_ties="all"' output to 'resolve_ties="first"'
blast <- blast %>%
    dplyr::group_by(qseqid) %>% 
    dplyr::mutate(row_n = dplyr::row_number()) %>%
    dplyr::top_n(1, row_n) %>% # Break ties by position
    dplyr::select(-row_n) %>%
    dplyr::ungroup()

### run R code

#filter tax table
tax <- tax %>% 
    seqateurs::unclassified_to_na(rownames=FALSE) %>%
    dplyr::mutate(lowest = seqateurs::lowest_classified(.)) # causes warning: ' argument is not an atomic vector; coercing'
    ### TODO: resolve above warning message

# make seqmap 
seqmap <- tibble::enframe(getSequences(seqtab), name = NULL, value="OTU") %>%
            mutate(name = paste0("SV", seq(length(getSequences(seqtab)))))
# get ASV sequences
seqs <- taxreturn::char2DNAbin(seqmap$OTU)
names(seqs) <- seqmap$name

# add sequences to blast output and rename columns
blast <- blast %>% 
    mutate(blastspp = paste0(Genus, " ", Species)) %>%
    dplyr::select(name = qseqid, acc, blastspp, pident, total_score, max_score, evalue, qcovs) %>%
    left_join(seqmap) %>%
    dplyr::select(-name)

# combine blast output and tax table
if ( nrow(blast) > 0 & nrow(tax) > 0 ){
                joint <- blast %>% 
                    dplyr::left_join(tax, by="OTU")
                } else {
                    NULL
                }

# save joint object
saveRDS(joint, paste0(fcid,"_",pcr_primers,"_joint.rds"))

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
            title = paste0(fcid, "  ", pcr_primers, " Top hit identity distribution"),
            subtitle = paste0("IDTAXA database:", idtaxa_db, " BLAST database:", ref_fasta),
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
            panel.grid = element_line(size = rel(0.5))
        ) 
} else { NULL }

# write plot
plot_filename <- paste0(fcid,"_",pcr_primers,"_taxonomic_assignment_summary.pdf")
if (!is.null(plot) ){
    pdf(file = plot_filename, width = 11, height = 8 , paper="a4r")
    print(plot)
    try(dev.off(), silent=TRUE)
} else {
    pdf(file = plot_filename, width = 11, height = 8 , paper="a4r")
    plot.new()
    text(x=.5, y=.5, "ERROR: No blast hits to reference fasta -- assignment plot not created") 
    try(dev.off(), silent=TRUE)
}

# stop(" *** stopped manually *** ") ##########################################