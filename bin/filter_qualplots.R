#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

reads_paths                 <- args$reads_paths
sample_primers              <- args$sample_primers
locus                       <- args$locus
primers                     <- args$primers
read_group                  <- args$read_group
seq_type                    <- args$seq_type
paired                      <- args$paired
file_suffix                 <- args$file_suffix

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "ggplot2",
    "patchwork",
    "seqateurs",
    "ShortRead",
    "stringr",
    "tidyr",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

### process variables 
# split reads_paths into individual file paths (or keep if single reads)
if ( paired == "true" ) {
    reads_paths_vec <- reads_paths %>% stringr::str_split_1(";")
    fwd_reads <- reads_paths_vec[1]
    rev_reads <- reads_paths_vec[2]
} else if ( paired == "false" ) {
   single_reads <- reads_paths
} else {
    stop ( "'paired' must be 'true' or 'false'!" )
}

truncLen <- NULL
quiet <- FALSE
n_reads <- 10000

### run R code

if ( paired == "true" & seq_type == "illumina" ){

    ### paired-end reads
    
    fastqFs <- normalizePath(fwd_reads)
    fastqRs <- normalizePath(rev_reads)

    if(!file.size(fastqFs) > 28) {
        message(paste0("Sample ", fwd_reads, " has no reads"))
        if ( file_suffix == "pre" ){
            empty_plot <- ggplot() + theme_void() + annotate(geom = "text", label = paste0("Sample '",sample_primers,"'\ncontains no reads pre-filtering"), x = 5, y = 5)
        } else if ( file_suffix == "post" ) {
            empty_plot <- ggplot() + theme_void() + annotate(geom = "text", label = paste0("Sample '",sample_primers,"'\ncontains no reads post-filtering"), x = 5, y = 5)
        } else {
            stop(paste0("'file_suffix' must be 'pre' or 'post'"))
        }
        gg.Fqual <- empty_plot
        gg.Rqual <- empty_plot
        gg.Fee <- empty_plot
        gg.Ree <- empty_plot

    } else {
        Fquals <- seqateurs::get_qual_stats(fastqFs, n=n_reads)
        Rquals <- seqateurs::get_qual_stats(fastqRs, n=n_reads)
            
        #Plot qualities
        gg.Fqual <- Fquals %>% 
            dplyr::select(Cycle, reads, starts_with("Q")) %>% 
            tidyr::pivot_longer(cols = starts_with("Q")) %>% 
            ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
            geom_line(data = Fquals, aes(y = QMean), color = "#66C2A5") + 
            geom_line(data = Fquals, aes(y = Q25), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
            geom_line(data = Fquals, aes(y = Q50), color = "#FC8D62", linewidth = 0.25) + 
            geom_line(data = Fquals, aes(y = Q75), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
            labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
            ggtitle(paste0(sample_primers," ",locus, " Forward Reads")) +
            theme(legend.position = "bottom", plot.title = element_text(size = 10))
        
        gg.Rqual <- Rquals %>% 
            dplyr::select(Cycle, reads, starts_with("Q")) %>% 
            tidyr::pivot_longer(cols = starts_with("Q")) %>% 
            ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
            geom_line(data = Rquals, aes(y = QMean), color = "#66C2A5") + 
            geom_line(data = Rquals, aes(y = Q25), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
            geom_line(data = Rquals, aes(y = Q50), color = "#FC8D62", linewidth = 0.25) + 
            geom_line(data = Rquals, aes(y = Q75), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
            labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
            ggtitle(paste0(sample_primers," ",locus, " Reverse Reads")) +
            scale_x_continuous(breaks=seq(0,300,25)) +
            theme(legend.position = "bottom", plot.title = element_text(size = 10))
        
        
        gg.Fee <- Fquals %>% 
            dplyr::select(Cycle, starts_with("EE")) %>% 
            tidyr::pivot_longer(cols = starts_with("EE")) %>% 
            dplyr::group_by(name) %>%  # Remove EE and replace with percentage at end? - i.e lower 10%
            dplyr::mutate(cumsumEE = cumsum(value)) %>%
            ggplot2::ggplot(aes(x = Cycle, y = log10(cumsumEE), colour = name)) + 
            geom_point(size = 1) + 
            geom_hline(yintercept = log10(1), color = "red") + 
            geom_hline(yintercept = log10(2), color = "red") + 
            geom_hline(yintercept = log10(3), color = "red") + 
            geom_hline(yintercept = log10(5), color = "red") + 
            geom_hline(yintercept = log10(7), color = "red") + 
            geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
            geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") + 
            labs(x = "Reads position", y = "Log10 Cumulative expected errors",
                colour = "Read quantiles") + 
            ggtitle(paste0(sample_primers," ",locus, " Forward Reads")) +
            scale_x_continuous(breaks=seq(0,300,25)) +
            theme(legend.position = "bottom", plot.title = element_text(size = 10))

        
        gg.Ree <- Rquals %>% 
            dplyr::select(Cycle, starts_with("EE")) %>% 
            tidyr::pivot_longer(cols = starts_with("EE")) %>% 
            dplyr::group_by(name) %>% 
            dplyr::mutate(cumsumEE = cumsum(value)) %>%
            ggplot2::ggplot(aes(x = Cycle, y = log10(cumsumEE), colour = name)) + 
            geom_point(size = 1) + 
            geom_hline(yintercept = log10(1), color = "red") + 
            geom_hline(yintercept = log10(2), color = "red") + 
            geom_hline(yintercept = log10(3), color = "red") + 
            geom_hline(yintercept = log10(5), color = "red") + 
            geom_hline(yintercept = log10(7), color = "red") + 
            geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
            geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") + 
            geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") + 
            labs(x = "Reads position", y = "Log10 Cumulative expected errors") +
            ggtitle(paste0(sample_primers," ",locus," Reverse Reads")) +
            scale_x_continuous(breaks=seq(0,300,25)) +
            theme(legend.position = "bottom", plot.title = element_text(size = 10))
        
        if(!is.null(truncLen)){
            gg.Fqual <- gg.Fqual +
            geom_vline(aes(xintercept=truncLen[1]), colour="blue") +
            annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
            gg.Fee <-gg.Fee +
            geom_vline(aes(xintercept=truncLen[1]), colour="blue")+
            annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
            gg.Rqual <- gg.Rqual +
            geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
            annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
            gg.Ree <- gg.Ree + 
            geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
            annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
        }
    }
    
    Qualplots <- (gg.Fqual + gg.Rqual) / (gg.Fee + gg.Ree)

    ggsave(paste0("qualplots_",file_suffix,"_",sample_primers,".pdf"), Qualplots, width = 11, height = 8)

} else if ( paired == "false" & seq_type == "nanopore" ) {

    ### single-end Nanopore reads
    fastqSs <- normalizePath(single_reads)

    if(!file.size(fastqSs) > 28) {
        message(paste0("Sample ", single_reads, " has no reads"))
        if ( file_suffix == "pre" ){
            empty_plot <- ggplot() + theme_void() + annotate(geom = "text", label = paste0("Sample '",sample_primers,"'\ncontains no reads pre-filtering"), x = 5, y = 5)
        } else if ( file_suffix == "post" ) {
            empty_plot <- ggplot() + theme_void() + annotate(geom = "text", label = paste0("Sample '",sample_primers,"'\ncontains no reads post-filtering"), x = 5, y = 5)
        } else {
            stop(paste0("'file_suffix' must be 'pre' or 'post'"))
        }
        gg.Squal <- empty_plot
        gg.See <- empty_plot
    }
    
    Squals <- seqateurs::get_qual_stats(fastqSs, n=n_reads)
        
    #Plot qualities
    gg.Squal <- Squals %>% 
        dplyr::select(Cycle, reads, starts_with("Q")) %>% 
        tidyr::pivot_longer(cols = starts_with("Q")) %>% 
        ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
        geom_line(data = Squals, aes(y = QMean), color = "#66C2A5") + 
        geom_line(data = Squals, aes(y = Q25), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
        geom_line(data = Squals, aes(y = Q50), color = "#FC8D62", linewidth = 0.25) + 
        geom_line(data = Squals, aes(y = Q75), color = "#FC8D62", linewidth = 0.25, linetype = "dashed") + 
        labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
        ggtitle(paste0(sample_primers," ",locus, " Single-end Reads")) +
        theme(legend.position = "bottom", plot.title = element_text(size = 10))
    
    gg.See <- Squals %>% 
        dplyr::select(Cycle, starts_with("EE")) %>% 
        tidyr::pivot_longer(cols = starts_with("EE")) %>% 
        dplyr::group_by(name) %>%  # Remove EE and replace with percentage at end? - i.e lower 10%
        dplyr::mutate(cumsumEE = cumsum(value)) %>%
        ggplot2::ggplot(aes(x = Cycle, y = log10(cumsumEE), colour = name)) + 
        geom_point(size = 1) + 
        geom_hline(yintercept = log10(1), color = "red") + 
        geom_hline(yintercept = log10(2), color = "red") + 
        geom_hline(yintercept = log10(3), color = "red") + 
        geom_hline(yintercept = log10(5), color = "red") + 
        geom_hline(yintercept = log10(7), color = "red") + 
        geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") + 
        geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
        geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") + 
        geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") + 
        geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") + 
        labs(x = "Reads position", y = "Log10 Cumulative expected errors",
            colour = "Read quantiles") + 
        ggtitle(paste0(sample_primers," ",locus, " Single-end Reads")) +
        # scale_x_continuous(breaks=seq(0,300,25)) +
        theme(legend.position = "bottom", plot.title = element_text(size = 10))

    # if(!is.null(truncLen)){
    #     gg.Fqual <- gg.Fqual +
    #     geom_vline(aes(xintercept=truncLen[1]), colour="blue") +
    #     annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    #     gg.Fee <-gg.Fee +
    #     geom_vline(aes(xintercept=truncLen[1]), colour="blue")+
    #     annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    #     gg.Rqual <- gg.Rqual +
    #     geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
    #     annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
    #     gg.Ree <- gg.Ree + 
    #     geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
    #     annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
    # }
    
    Qualplots <- gg.Squal / gg.See

    ggsave(paste0("qualplots_",file_suffix,"_",sample_primers,".pdf"), Qualplots, width = 11, height = 8)

} else {
    stop ( "Currently unsupported sequencing type -- check samplesheet" )
}

# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})