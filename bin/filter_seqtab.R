#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dada2",
    "dplyr",
    "ggplot2",
    "patchwork",
    "readr",
    "stringr",
    "taxreturn",
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
    "asv_min_length",
    "asv_max_length",
    "phmm",
    "coding",
    "genetic_code",
    "for_primer_seq",
    "rev_primer_seq"
)
lapply(nf_vars, nf_var_check)

## extract variables and deal with repetition
asv_min_length <-   parse_nf_var_repeat(asv_min_length) %>% as.numeric
asv_max_length <-   parse_nf_var_repeat(asv_max_length) %>% as.numeric
phmm <-             parse_nf_var_repeat(phmm)
coding <-           parse_nf_var_repeat(coding) %>% as.logical
genetic_code <-     parse_nf_var_repeat(genetic_code)
for_primer_seq <-   parse_nf_var_repeat(for_primer_seq)
rev_primer_seq <-   parse_nf_var_repeat(rev_primer_seq)

## check and define variables 
if(is.na(asv_min_length))   {asv_min_length <- NULL}
if(is.na(asv_max_length))   {asv_max_length <- NULL}
if(is.na(phmm))             {phmm <- NULL}
if(is.na(for_primer_seq))   {for_primer_seq <- NULL}
if(is.na(rev_primer_seq))   {rev_primer_seq <- NULL}

check_frame <- coding
quiet <- FALSE # switch quiet off for now
multithread <- FALSE # multithreading switched off for now
### TODO: Implement multithreading

# set primers vector
primers <- c(for_primer_seq, rev_primer_seq)

# read seqtab from file 
if(is.matrix(seqtab) | is.data.frame(seqtab)){
    if(!quiet){message("Input is a matrix or data frame")}
} else if (is.character(seqtab) & stringr::str_detect(seqtab, ".rds")){
    seqtab <- readRDS(seqtab)
} else {
    stop("seqtab must be a matrix/data frame or .rds file")
}

### run R code

reads_starting <- rowSums(seqtab)

## Load in profile hidden markov model if provided
if(is.character(phmm) && stringr::str_detect(phmm, ".rds")){
    phmm_model <- readRDS(phmm)
    message("Loaded PHMM from file")
} else if (is(phmm, "PHMM")){
    phmm_model <- phmm
    message("Loaded PHMM from R object")
} else {
    phmm_model <- NULL
    message("Running analysis with no PHMM")
}

## subset PHMM if primers were provided
if (is(phmm_model, "PHMM") && !is.null(primers)){
    # Check that one of the two primers can bind
    Fbind <- taxreturn::get_binding_position(primers[1], model = phmm_model, tryRC = TRUE, min_score = 10)
    Rbind <- taxreturn::get_binding_position(primers[2], model = phmm_model, tryRC = TRUE, min_score = 10)
    if(!is.na(Fbind$start) & !is.na(Rbind$start)){
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else if(!is.na(Fbind$start) & is.na(Rbind$start)){
        # Reverse primer not found - Try with subsets
        for(r in seq(1, nchar(primers[2])-10, 1)){ #Minimum length of 10 as this has to match minscore
            Rbind <- get_binding_position(
                str_remove(primers[2], paste0("^.{1,",r,"}")), 
                model = phmm_model, 
                tryRC = TRUE, 
                min_score = 10
                )
            if (!is.na(Rbind$start)) {
                primers[2] <- stringr::str_remove(primers[2], paste0("^.{1,",r,"}"))
                break
            }
        }
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else  if(is.na(Fbind$start) & !is.na(Rbind$start)){
        # Forward primer not found - Try with subsets
        for(r in seq(1, nchar(primers[1])-10, 1)){ #Minimum length of 10 as this has to match minscore
        Rbind <- taxreturn::get_binding_position(stringr::str_remove(primers[1], paste0("^.{1,",r,"}")), model = phmm_model, tryRC = TRUE, min_score = 10)
        if (!is.na(Rbind$start)) {
            primers[1] <- stringr::str_remove(primers[1], paste0("^.{1,",r,"}"))
            break
        }
        }
        phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    }
}

## Remove chimeras  
seqtab_nochim <- dada2::removeBimeraDenovo(
    seqtab, method="consensus", multithread=multithread, verbose=!quiet
    )
seqs_rem <- length(colnames(seqtab_nochim))/length(colnames(seqtab))
abund_rem <- sum(seqtab_nochim)/sum(seqtab)
message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% initial abundance remaining after chimera removal"))
reads_chimerafilt <- rowSums(seqtab_nochim)

## subset ASVs to those of expected size
if(any(!is.null(c(asv_min_length, asv_max_length)), na.rm = TRUE) & any(reads_chimerafilt > 0)){
    if(!is.null(asv_min_length) & !is.null(asv_max_length)){
        seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) < asv_max_length & nchar(colnames(seqtab_nochim)) > asv_min_length, drop = FALSE]
    } else if(is.null(asv_min_length) & !is.null(asv_max_length)){
        seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) < asv_max_length, drop = FALSE]
    } else if(!is.null(asv_min_length) & is.null(asv_max_length)){
        seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) > asv_min_length, drop = FALSE]
    }
    if(!quiet){
        seqs_rem <- length(colnames(seqtab_cut))/length(colnames(seqtab))
        abund_rem <- sum(seqtab_cut)/sum(seqtab)
        message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% initial abundance remaining after length filtering"))
    }
    reads_lengthfilt <- rowSums(seqtab_cut)
} else {
    seqtab_cut <- seqtab_nochim
    reads_lengthfilt <- rep(0,nrow(seqtab)) 
    names(reads_lengthfilt) <- rownames(seqtab)
}

## Align ASVs against phmm
if (is(phmm_model, "PHMM") & any(reads_lengthfilt > 0)){
    seqs <- Biostrings::DNAStringSet(colnames(seqtab_cut))
    names(seqs) <- colnames(seqtab_cut)
    phmm_filt <- taxreturn::map_to_model(
        seqs, model = phmm_model, min_score = 100, min_length = 100,
        shave = FALSE, check_frame = check_frame, kmer_threshold = 0.5, k=5, extra = "fill")
    seqtab_phmm <- seqtab_cut[,colnames(seqtab_cut) %in% names(phmm_filt), drop = FALSE]
    if(!quiet){
        seqs_rem <- length(colnames(seqtab_phmm))/length(colnames(seqtab))
        abund_rem <- sum(seqtab_phmm)/sum(seqtab)
        message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% of initial abundance remaining after PHMM filtering"))
    }
    reads_phmmfilt <- rowSums(seqtab_phmm)
} else {
    seqtab_phmm <- seqtab_cut
    reads_phmmfilt <- rep(0,nrow(seqtab)) 
    names(reads_phmmfilt) <- rownames(seqtab)
}

print(reads_phmmfilt)
print(check_frame)

# stop(" *** stopped manually *** ") ##########################################

## Filter sequences containing stop codons
if(check_frame & any(reads_phmmfilt > 0)){
    seqs <- Biostrings::DNAStringSet(colnames(seqtab_phmm))
    names(seqs) <- colnames(seqtab_phmm)
    codon_filt <- taxreturn::codon_filter(seqs, genetic_code = genetic_code) 
    seqtab_final <- seqtab_phmm[,colnames(seqtab_phmm) %in% names(codon_filt), drop = FALSE]
    if(!quiet){
        seqs_rem <- length(colnames(seqtab_final))/length(colnames(seqtab))
        abund_rem <- sum(seqtab_final)/sum(seqtab)
        message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% of initial abundance remaining after checking reading frame"))
    }
    reads_framefilt <- rowSums(seqtab_final)
} else {
    seqtab_final <- seqtab_phmm
    reads_framefilt <- rep(0,nrow(seqtab))
    names(reads_framefilt) <- rownames(seqtab)
}

# change sequence names to contain primer names (preceded by single underscore), so merging seqtabs across loci with sample sheets works as intended
seqtab_final_renamed <- seqtab_final %>% 
    tibble::as_tibble(rownames = "sample_id") %>%
    dplyr::mutate(sample_id = paste0(sample_id,"_",pcr_primers)) %>%
    tibble::column_to_rownames("sample_id") %>% 
    as.matrix()

reads_final <- rowSums(seqtab_final_renamed)

## Output a cleanup summary
cleanup <- seqtab %>%
    as.data.frame() %>%
    tidyr::pivot_longer( tidyselect::everything(),
                    names_to = "OTU",
                    values_to = "Abundance") %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise(Abundance = sum(Abundance)) %>%
    dplyr::mutate(length  = nchar(OTU)) %>%
    dplyr::mutate(type = dplyr::case_when(
        !OTU %in% dada2::getSequences(seqtab_nochim) ~ "Chimera",
        !OTU %in% dada2::getSequences(seqtab_cut) ~ "Incorrect size",
        !OTU %in% dada2::getSequences(seqtab_phmm) ~ "PHMM",
        !OTU %in% dada2::getSequences(seqtab_final) ~ "Stop codons",
        TRUE ~ "Retained"
    )) %>%
    dplyr::mutate(concat = str_detect(OTU, "NNNNNNNNNN")) %>%
    dplyr::mutate(type = case_when(
        concat ~ paste0("Unmerged-", type),
        TRUE ~ type
    )) %>%
    dplyr::mutate(pcr_primers = pcr_primers) %>%
    dplyr::select(-concat)

## Output length distribution plots

# colours of columns
cols <- c(`Chimera` = "#9e0142",
        `Unmerged-Chimera` = "#d53e4f",
        `Incorrect size` = "#f46d43",
        `Unmerged-Incorrect size` = "#fdae61",
        `PHMM` = "#fee08b",
        `Unmerged-PHMM` = "#e6f598",
        `Stop codons` = "#abdda4",
        `Unmerged-Stop codons` = "#66c2a5",
        `Retained` = "#3288bd",
        `Unmerged-Retained` = "#5e4fa2") 

gg.abundance <- 
    ggplot2::ggplot(cleanup, aes(x=length, y=log10(Abundance), fill=type)) +
    geom_col() + 
    # scale_y_continuous(trans = pseudo_log_trans(sigma = 1)) + 
    scale_x_continuous(limits=c(min(cleanup$length)-10, max(cleanup$length)+10))+
    theme_bw()+
    scale_fill_manual(values = cols)+
    theme(
        strip.background = element_rect(colour = "black", fill = "lightgray"),
        strip.text = element_text(size=9, family = ""),
        axis.text.x =element_text(angle=45, hjust=1, vjust=1),
        plot.background = element_blank(),
        text = element_text(size=9, family = ""),
        axis.text = element_text(size=8, family = ""),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_line(size = rel(0.5)),
    ) +
    labs(
        title=pcr_primers,
        subtitle = "Abundance of sequences",
        x = "ASV length",
        y = "log10 ASV abundance",
        fill = "ASV type"
        )

gg.unique <- 
    ggplot2::ggplot(cleanup, aes(x=length, fill=type))+
    geom_histogram(binwidth = 1) + 
    scale_x_continuous(limits=c(min(cleanup$length)-10, max(cleanup$length)+10))+
    theme_bw()+
    scale_fill_manual(values = cols)+
    theme(
        strip.background = element_rect(colour = "black", fill = "lightgray"),
        strip.text = element_text(size=9, family = ""),
        axis.text.x =element_text(angle=45, hjust=1, vjust=1),
        plot.background = element_blank(),
        text = element_text(size=9, family = ""),
        axis.text = element_text(size=8, family = ""),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid = element_line(size = rel(0.5)),
    ) +
    labs(
        title=pcr_primers,
        subtitle = "Number of unique sequences",
        x = "ASV length",
        y = "Number of unique sequences",
        fill = "ASV type"
        )

## Create combined plot
out_plot <- gg.abundance / gg.unique


## outputs
# plot
ggsave(paste0(fcid,"_",pcr_primers,"_ASV_cleanup_summary.pdf"), out_plot, width = 11, height = 8)
# summary table
dplyr::bind_rows(unique(cleanup)) %>%
    readr::write_csv(paste0(fcid,"_",pcr_primers,"_ASV_cleanup_summary.csv"))
# filtered seqtab
saveRDS(seqtab_final_renamed, paste0(fcid,"_",pcr_primers,"_seqtab.cleaned.rds"))

# for read-tracking
reads_out <- tibble::tibble(
    sample_id = rownames(seqtab) %>% stringr::str_remove(pattern="_S[0-9]+_R[1-2]_.*$"),
    filter_chimera = reads_chimerafilt,
    filter_length = reads_lengthfilt,
    filter_phmm = reads_phmmfilt,
    filter_frame = reads_framefilt,
    filter_seqtab = reads_final
    ) %>%
    tidyr::pivot_longer(
        cols = filter_chimera:filter_seqtab,
        names_to = "stage",
        values_to = "pairs"
        ) %>% 
    dplyr::mutate(
        fcid = fcid, 
        pcr_primers = pcr_primers
    ) %>%
    dplyr::select(stage, sample_id, fcid, pcr_primers, pairs)

readr::write_csv(reads_out, paste0("filter_seqtab_",fcid,"_",pcr_primers,"_readsout.csv"))

# stop(" *** stopped manually *** ") ##########################################

