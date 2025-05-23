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
    "asv_min_length",
    "asv_max_length",
    "phmm",
    "coding",
    "genetic_code",
    "for_primer_seq",
    "rev_primer_seq",
    "seqtab_tibble",
    "fasta"
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

seqtab_tibble <- readr::read_csv(seqtab_tibble)

seqs_fasta <- Biostrings::readDNAStringSet(fasta)

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


### run R code

# convert DSS to tibble of name and seq
seqs_names <- 
    seqs_fasta %>% 
    as.character() %>% 
    tibble::enframe(name = "seq_name", value = "sequence")

# make a dada2-style seqtab from inputs
seqtab_matrix <- 
    seqtab_tibble %>%
    # add seq string
    dplyr::left_join(., seqs_names, by = "seq_name") %>%
    # remove seq_name
    dplyr::select(-seq_name) %>%
    # pivot longer
    tidyr::pivot_longer(cols = !sequence, names_to = "sample_id", values_to = "abundance") %>%
    dplyr::mutate(abundance = as.integer(abundance)) %>%
    # pivot wider 
    tidyr::pivot_wider(names_from = sequence, values_from = abundance) %>%
    # sample_id as rownames
    tibble::column_to_rownames(var = "sample_id") %>%
    # convert to matrix
    as.matrix() 


## Remove chimeras  
seqtab_nochim <- 
    dada2::removeBimeraDenovo(
		seqtab_matrix, 
		method = "consensus", 
		multithread = multithread, 
		verbose = !quiet
    )

# make tibble of name, sequence and whether it passed the chimera filter
seqs_chimera_filter <- 
	seqs_names %>%
	dplyr::mutate(
		chimera_filter = sequence %in% colnames(seqtab_nochim)
	)


## remove based on length
# use seqs_names for this as easier
seqs_length_filter <- 
    seqs_names %>%
    dplyr::mutate(
        length = stringr::str_count(sequence, "\\w"),
        pass_min = length > asv_min_length,
        pass_max = length < asv_max_length,
        length_filter = pass_min & pass_max
    ) %>%
    # make tibble of name, sequence, and whether it passed the length filter
    dplyr::select(-length, -pass_min, -pass_max)

## PHMM filtering
if (is(phmm_model, "PHMM")){
    # remove sequences that don't align to PHMM
    phmm_filt_seqtab <- 
        taxreturn::map_to_model(
            seqs_fasta, 
            model = phmm_model, 
            min_score = 100, 
            min_length = 100,
            shave = FALSE, 
            check_frame = check_frame, 
            kmer_threshold = 0.5, 
            k = 5, 
            extra = "fill"
        )
    
    # make tibble of name, sequence and whether it passed the PHMM filter
    seqs_phmm_filter <- 
        seqs_names %>%
        dplyr::mutate(
            phmm_filter = seq_name %in% names(phmm_filt_seqtab)
        )

} else {
    # all seqs pass filter
    seqs_phmm_filter <- 
        seqs_names %>%
        dplyr::mutate(
            phmm_filter = TRUE
        )
}

## frame checking
if(check_frame){
    # remove sequences with wrong frame  
    codon_filt_dss <- 
        taxreturn::codon_filter(
            seqs_fasta, 
            genetic_code = genetic_code
        ) 
    
    # make tibble of name, sequence and whether it passed the PHMM filter
    seqs_frame_filter <- 
        seqs_names %>%
        dplyr::mutate(
            frame_filter = seq_name %in% names(codon_filt_dss)
        )
        
    } else {
    # all seqs pass filter
    seqs_frame_filter <- 
        seqs_names %>%
        dplyr::mutate(
            frame_filter = TRUE
        )
}

## combine all the filters together and add to seqtab_tibble
seqs_names_filters <- 
    seqs_names %>%
    dplyr::left_join(., seqs_chimera_filter, by = dplyr::join_by(seq_name, sequence)) %>%
    dplyr::left_join(., seqs_length_filter, by = dplyr::join_by(seq_name, sequence)) %>%
    dplyr::left_join(., seqs_phmm_filter, by = dplyr::join_by(seq_name, sequence)) %>%
    dplyr::left_join(., seqs_frame_filter, by = dplyr::join_by(seq_name, sequence)) %>%
    # join to seqtab_tibble
    dplyr::left_join(seqtab_tibble, ., by = dplyr::join_by(seq_name))

# save seqtab_tibble
readr::write_csv(seqs_names_filters, paste0(fcid, "_", pcr_primers, "_seqtab_filtered.csv"))

# filter info tibble

filter_info <- 
    seqs_names_filters %>%
    tidyr::pivot_longer(
        cols = !c(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter), 
        names_to = "sample_id",
        values_to = "abundance"
    ) %>%
    # summarise down to total abundance across all samples (one row per ASV)
    dplyr::group_by(seq_name) %>%
    dplyr::mutate(abundance = sum(abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sample_id) %>%
    dplyr::distinct() %>%
    # add extra info
    dplyr::mutate(
        length = nchar(sequence),
        pcr_primers = pcr_primers,
        concat = stringr::str_detect(sequence, rep(x = "N", times = 10) %>% stringr::str_flatten()),
        dplyr::across(
            chimera_filter:frame_filter, 
            ~ dplyr::case_when(
                . == TRUE ~ "pass", 
                . == FALSE ~ "fail"
            )
        ),
        dplyr::across(
            chimera_filter:frame_filter, 
            ~ forcats::fct_relevel(., "pass", "fail")
        )
    )

# save filter summary
readr::write_csv(filter_info, paste0(fcid, "_", pcr_primers, "_ASV_cleanup.csv"))

## plots

# Output length distribution plots

# function to make generic plot
ab_plot_fun <- 
    function(
        filter_column,
        filter_name,
        filter_colour,
        dataset = filter_info
    ){
        filter_column <- enquo(filter_column)
        # summarise filter
        dataset %>%
            dplyr::group_by(!!filter_column, length, concat) %>%
            dplyr::summarise(abundance = sum(abundance)) %>%
            dplyr::ungroup() %>%
            # plot
            ggplot2::ggplot(., aes(x = length, y = abundance, colour = !!filter_column, shape = concat)) +
            geom_point(size = 3, alpha = 0.7) + 
            scale_y_continuous(
                trans = scales::pseudo_log_trans(sigma = 1), 
                limits = c(0, dataset$abundance %>% sum),
                breaks = c(0, 10^seq(1,100,1))
            ) +
            scale_x_continuous(limits=c(min(dataset$length)-10, max(dataset$length)+10)) +
            theme_bw() +
            scale_colour_manual(values = c("pass" = "grey80", "fail" = filter_colour)) +
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
                subtitle = filter_name,
                x = "ASV length",
                y = "ASV abundance",
                fill = filter_name
            )
    }

# make plots
gg.chimera.ab <- 
    ab_plot_fun(
        chimera_filter,
        "Chimera filter",
        "#9e0142"
    )

gg.length.ab <- 
    ab_plot_fun(
        length_filter,
        "Length filter",
        "#f46d43"
    )
    
gg.phmm.ab <-
    ab_plot_fun(
        phmm_filter,
        "PHMM filter",
        "#fee08b"
    )
    
gg.frame.ab <- 
    ab_plot_fun(
        frame_filter,
        "Frame/stop filter",
        "#abdda4"
    )  

# combine abundance plots

gg.abundance <-
    gg.chimera.ab + gg.length.ab + gg.phmm.ab + gg.frame.ab +
    patchwork::plot_layout(ncol = 1, guides = "collect") +
    patchwork::plot_annotation(title = "ASV abundance", subtitle = paste0("Flowcell: ", fcid, "\nPCR primers: ", pcr_primers))

ggsave(paste0(fcid, "_", pcr_primers,"_asv_abundance.pdf"), gg.abundance, width = 8, height = 12)


## ASV count plots

# function to make generic plot
n_plot_fun <- 
    function(
        filter_column,
        filter_name,
        filter_colour,
        dataset = filter_info
    ){
        filter_column <- enquo(filter_column)
        # summarise filter
        dataset %>%
            dplyr::group_by(!!filter_column, length, concat) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::ungroup() %>%
            # plot
            ggplot2::ggplot(., aes(x = length, y = n, colour = !!filter_column, shape = concat)) +
            geom_point(size = 3) + 
            scale_y_continuous(
                limits = c(0, dataset$seq_name %>% unique() %>% length())
            ) +
            scale_x_continuous(limits=c(min(dataset$length)-10, max(dataset$length)+10)) +
            theme_bw() +
            scale_colour_manual(values = c("pass" = "grey80", "fail" = filter_colour)) +
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
                subtitle = filter_name,
                x = "ASV length",
                y = "ASV count",
                fill = filter_name
            )
    }

# make plots
gg.chimera.n <- 
    n_plot_fun(
        chimera_filter,
        "Chimera filter",
        "#9e0142"
    )

gg.length.n <- 
    n_plot_fun(
        length_filter,
        "Length filter",
        "#f46d43"
    )
    
gg.phmm.n <-
    n_plot_fun(
        phmm_filter,
        "PHMM filter",
        "#fee08b"
    )
    
gg.frame.n <- 
    n_plot_fun(
        frame_filter,
        "Frame/stop filter",
        "#abdda4"
    )  

# combine count plots

gg.n_asvs <-
    gg.chimera.n + gg.length.n + gg.phmm.n + gg.frame.n +
    patchwork::plot_layout(ncol = 1, guides = "collect") +
    patchwork::plot_annotation(title = "Unique ASVs", subtitle = paste0("Flowcell: ", fcid, "\nPCR primers: ", pcr_primers))

ggsave(paste0(fcid, "_", pcr_primers,"_asv_count.pdf"), gg.n_asvs, width = 8, height = 12)


## for read-tracking

# how many reads pass through each filter for each sample?
reads_out <- 
    seqs_names_filters %>% 
    # rename filters to match format required for read tracking
    dplyr::rename(
        "filter_chimera" = chimera_filter,
        "filter_length" = length_filter,
        "filter_phmm" = phmm_filter,
        "filter_frame" = frame_filter
    ) %>%
    # add "filter_seqtab" which is passed when all filters are passed
    dplyr::mutate(filter_seqtab = filter_chimera & filter_length & filter_phmm & filter_frame) %>%
    tidyr::pivot_longer(
        cols = !c(seq_name, sequence, filter_chimera, filter_length, filter_phmm, filter_frame, filter_seqtab),
        names_to = "sample_id",
        values_to = "abundance"
    ) %>%
    dplyr::select(-seq_name, -sequence) %>%
    dplyr::mutate(fcid = fcid, pcr_primers = pcr_primers) %>%
    tidyr::pivot_longer(
        cols = tidyselect::starts_with("filter_"), 
        names_to = "stage",
        values_to = "pass"
    ) %>%
    dplyr::group_by(sample_id, fcid, pcr_primers, stage, pass) %>%
    dplyr::summarise(pairs = sum(abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(pass == TRUE) %>%
    dplyr::select(-pass) %>%
    dplyr::select(stage, sample_id, fcid, pcr_primers, pairs)

readr::write_csv(reads_out, paste0("filter_seqtab_",fcid,"_",pcr_primers,"_readsout.csv"))
