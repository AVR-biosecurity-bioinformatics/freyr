#!/usr/bin/env Rscript

# check variables defined


### run R code

plot_read_quals2 <- function(sample_id, fwd_reads, rev_reads, fcid, target_gene, pcr_primers, truncLen = NULL, quiet=FALSE, n = 10000){
    fastqFs <- fwd_reads
    fastqRs <- rev_reads
    
    if(length(fastqFs) == 0 ){
        message(paste0("Sample ", fwd_reads, " Has no reads"))
        return(NULL)
    }
    if(!file.size(fastqFs) > 28) {
        message(paste0("Sample ", fwd_reads, " Has no reads"))
        return(NULL)
    }
    
    Fquals <- get_qual_stats(fastqFs, n=n)
    Rquals <- get_qual_stats(fastqRs, n=n)
    
    #Plot qualities
    gg.Fqual <- Fquals %>% 
        dplyr::select(Cycle, reads, starts_with("Q")) %>% 
        tidyr::pivot_longer(cols = starts_with("Q")) %>% 
        ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
        geom_line(data = Fquals, aes(y = QMean), color = "#66C2A5") + 
        geom_line(data = Fquals, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        geom_line(data = Fquals, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
        geom_line(data = Fquals, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
        ggtitle(paste0(sample_id," ",target_gene," ",pcr_primers, " Forward Reads")) +
        scale_x_continuous(breaks=seq(0,300,25))
    
    gg.Rqual <- Rquals %>% 
        dplyr::select(Cycle, reads, starts_with("Q")) %>% 
        tidyr::pivot_longer(cols = starts_with("Q")) %>% 
        ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
        geom_line(data = Rquals, aes(y = QMean), color = "#66C2A5") + 
        geom_line(data = Rquals, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        geom_line(data = Rquals, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
        geom_line(data = Rquals, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
        ggtitle(paste0(sample_id," ",target_gene," ",pcr_primers, " Reverse Reads")) +
        scale_x_continuous(breaks=seq(0,300,25)) +
        theme(legend.position = "bottom")
    
    
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
        ggtitle(paste0(sample_id," ",target_gene," ",pcr_primers, " Forward Reads")) +
        scale_x_continuous(breaks=seq(0,300,25)) +
        theme(legend.position = "bottom")

    
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
        ggtitle(paste0(sample_id," ",target_gene," ",pcr_primers, " Reverse Reads")) +
        scale_x_continuous(breaks=seq(0,300,25)) +
        theme(legend.position = "bottom")
    
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
    
    Qualplots <- (gg.Fqual + gg.Rqual) / (gg.Fee + gg.Ree)

    ggsave(Qualplots, paste0(sample_id,"_",target_gene,"_",pcr_primers,"_qualplots.pdf"))
    
    return(Qualplots)
}




plot_read_quals2(
    sample_id =     sample_id, 
    fwd_reads =     fwd_reads, 
    rev_reads =     rev_reads, 
    fcid =          fcid,
    target_gene =   target_gene,
    pcr_primers =   pcr_primers, 
    truncLen = NULL, 
    quiet=FALSE, 
    n = 1000
)
