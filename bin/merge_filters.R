#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

primers                     <- args$primers
filter_tibble_list          <- args$filter_tibble_list
seqtab_tibble_list          <- args$seqtab_tibble_list
fasta_list                  <- args$fasta_list
samplesheet_split           <- args$samplesheet_split

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "ggplot2",
    "patchwork",
    "readr",
    "stringr",
    "tibble",
    NULL
)
suppressPackageStartupMessages(invisible(lapply(process_packages, library, character.only = TRUE, warn.conflicts = FALSE)))

## check and define variables

filter_list <- 
    filter_tibble_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) # read in tibbles and store as list of tibbles

seqtab_list <- 
    seqtab_tibble_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., readr::read_csv) # read in seqtabs and store as list of tibbles

fasta_list <- 
    fasta_list %>%
    stringr::str_split_1(., pattern = " ") %>% # split string of filenames into vector
    lapply(., Biostrings::readDNAStringSet)

samplesheet_tibble <- readr::read_csv(samplesheet_split)

### run R code

# combine filters into one tibble (one row per sequence)
filters_combined <- 
    filter_list %>%
    lapply(
        ., 
        function(x){
            tidyr::pivot_longer(
                x, 
                cols = tidyselect::ends_with("_filter"), 
                names_to = "filter", 
                values_to = "status"
            )
        }
    ) %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_wider(names_from = filter, values_from = status)

# export combined filters
readr::write_csv(filters_combined, paste0(primers,"_filters.csv"))

# combine seqtabs together (wide format)
seqtab_combined <- 
    seqtab_list %>%
    # pivot each tibble longer
    lapply(
        .,
        function(x){ # per tibble
            x %>%
            tidyr::pivot_longer(
                cols = !seq_name,
                names_to = "sample_primers",
                values_to = "abundance"
            )
        }
    ) %>%
    # bind tibbles together now columns all match
    dplyr::bind_rows() %>%
    # pivot wider, filling missing abundance with 0
    tidyr::pivot_wider(
        names_from = sample_primers,
        values_from = abundance, 
        values_fill = 0
    )

rm(seqtab_list)
gc()

# combine fasta files together, removing redundancy, and convert to table
seqs_names <-   
    fasta_list %>%
    lapply(., as.character) %>%
    unlist(use.names = T) %>% 
    tibble::enframe(name = "seq_name", value = "sequence") %>%
    dplyr::group_by(seq_name, sequence) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 

rm(fasta_list)
gc()

# add sequence to seqtab and export
seqtab_combined %>%
    dplyr::left_join(., seqs_names, by = "seq_name") %>%
    dplyr::relocate(seq_name, sequence) %>%
    readr::write_csv(., paste0(primers, "_seqtab_combined.csv"))

# export combined .fasta of all sequences across flowcells 
seqs_names %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet() %>%
    write_fasta(., paste0(primers,"_combseqs.fasta"))

# get tibble of read_group and sample_primers
sample_read_group <- 
    samplesheet_tibble %>%
    dplyr::select(sample_primers, read_group)

## make new read tracking tibble from filters and seqtab
filter_tracking <- 
    seqs_names %>%
    dplyr::left_join(., filters_combined, by = dplyr::join_by(seq_name, sequence)) %>%
    # make combined filter that is passed when all filters are passed
    dplyr::mutate(combined_filter = chimera_filter & length_filter & phmm_filter & frame_filter) %>%
    # join seqtab
    dplyr::left_join(., seqtab_combined, by = dplyr::join_by(seq_name)) %>%
    # pivot longer
    tidyr::pivot_longer(
        cols = !c(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter, combined_filter), 
        names_to = "sample_primers",
        values_to = "abundance"
    ) %>%
    # add read_group 
    dplyr::left_join(
        ., 
        sample_read_group, 
        by = dplyr::join_by("sample_primers")
    )
  
# make list of tibbles for filter plotting (and export), one tibble per flowcell 
filter_info <- 
    filter_tracking %>%
    # summarise down to total abundance across all samples (one row per ASV)
    dplyr::group_by(seq_name, read_group) %>%
    dplyr::mutate(abundance = sum(abundance)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-sample_primers) %>%
    dplyr::distinct() %>%
    # add extra info
    dplyr::mutate(
        length = nchar(sequence),
        primers = primers,
        concat = stringr::str_detect(sequence, rep(x = "N", times = 10) %>% stringr::str_flatten()),
        dplyr::across(
            c(chimera_filter, length_filter, phmm_filter, frame_filter, combined_filter), 
            ~ dplyr::case_when(
                . == TRUE ~ "pass", 
                . == FALSE ~ "fail"
            )
        ),
        dplyr::across(
            c(chimera_filter, length_filter, phmm_filter, frame_filter, combined_filter), 
            ~ forcats::fct_relevel(., "pass", "fail")
        )
    ) %>%
    # remove rows with abundance = 0
    dplyr::filter(abundance > 0) %>%
    # enforce column order
    dplyr::select(seq_name, sequence, chimera_filter, length_filter, phmm_filter, frame_filter, combined_filter, abundance, length, read_group, primers, concat) %>%
    # arrange by abundance
    dplyr::arrange(desc(abundance)) %>%
    # split into list, one tibble per read_group
    split(., .$read_group)

# save filter summary files, one per flowcell
filter_info %>%
    lapply(
        ., 
        function(x){
            readr::write_csv(x, paste0("asv_cleanup_",unique(x$read_group), "_", unique(x$primers), ".csv"))
        }
    )

gc()

## export read tracking group tibble per flowcell (how many reads pass through each filter for each sample?)
filter_tracking %>%
    # split into list by read_group
    split(., .$read_group) %>%
    lapply(
        ., 
        function(x){
            x %>%
            tidyr::pivot_wider(
                names_from = sample_primers, 
                values_from = abundance
            ) %>%
            dplyr::rename(
                "filter_chimera" = chimera_filter,
                "filter_length" = length_filter,
                "filter_phmm" = phmm_filter,
                "filter_frame" = frame_filter,
                "filter_combined" = combined_filter
            ) %>%
            tidyr::pivot_longer(
                cols = !c(seq_name, sequence, filter_chimera, filter_length, filter_phmm, filter_frame, filter_combined, read_group),
                names_to = "sample_primers",
                values_to = "abundance"
            ) %>%
            dplyr::select(-seq_name, -sequence) %>%
            dplyr::mutate(primers = primers) %>%
            tidyr::pivot_longer(
                cols = tidyselect::starts_with("filter_"), 
                names_to = "stage",
                values_to = "pass"
            ) %>%
            dplyr::group_by(sample_primers, read_group, primers, stage, pass) %>%
            dplyr::summarise(pairs = sum(abundance)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(pass == TRUE) %>%
            dplyr::select(-pass) %>%
            dplyr::select(stage, sample_primers, read_group, primers, pairs) %>%
            readr::write_csv(., paste0("filter_merged_",unique(.$read_group),"_", primers,"_readsout.csv"))
        }
    )

gc()

### generate filtering plots

# function to make abundance plot
ab_plot_fun <- 
    function(
        filter_column,
        filter_name,
        filter_colour,
        dataset
    ){
        filter_column <- enquo(filter_column)
        # summarise filter
        dataset %>%
            # dplyr::group_by(!!filter_column, length, concat) %>%
            # dplyr::summarise(abundance = sum(abundance)) %>%
            # dplyr::ungroup() %>%
            # plot
            ggplot2::ggplot(., aes(x = length, y = abundance, colour = !!filter_column, shape = concat)) +
            geom_point(size = 2.5, alpha = 0.7) + 
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

# function to make count plot
n_plot_fun <- 
    function(
        filter_column,
        filter_name,
        filter_colour,
        dataset
    ){
        filter_column <- enquo(filter_column)
        # summarise filter
        dataset %>%
            dplyr::group_by(!!filter_column, length, concat) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::ungroup() %>%
            # plot
            ggplot2::ggplot(., aes(x = length, y = n, colour = !!filter_column, shape = concat)) +
            geom_point(size = 2.5, alpha = 0.7) + 
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

# generate and write one combined abundance plot per flowcell
filter_info %>%
    lapply(
        ., 
        function(x){
        # make plots
        gg.combined.ab <- 
            ab_plot_fun(
                combined_filter,
                "Combined filters",
                "black",
                x
            )
        
        gg.chimera.ab <- 
            ab_plot_fun(
                chimera_filter,
                "Chimera filter",
                "#9e0142",
                x
            )
        
        gg.length.ab <- 
            ab_plot_fun(
                length_filter,
                "Length filter",
                "#f46d43",
                x
            )
            
        gg.phmm.ab <-
            ab_plot_fun(
                phmm_filter,
                "PHMM filter",
                "#fee08b",
                x
            )
            
        gg.frame.ab <- 
            ab_plot_fun(
                frame_filter,
                "Frame/stop filter",
                "#abdda4",
                x
            )  
        
        # combine abundance plots
        
        gg.abundance <-
            gg.combined.ab + gg.chimera.ab + gg.length.ab + gg.phmm.ab + gg.frame.ab +
            patchwork::plot_layout(ncol = 1, guides = "collect") +
            patchwork::plot_annotation(
                title = "ASV abundance", 
                subtitle = paste0("Flowcell: ", unique(x$read_group), "\nPCR primers: ", unique(x$primers))
            )
        
        ggsave(paste0("asv_abundance_",unique(x$read_group), "_", unique(x$primers),".pdf"), gg.abundance, width = 8, height = 12)

        }
    )

gc()

# generate and write one combined count plot per flowcell
filter_info %>%
    lapply(
        ., 
        function(x){
        # make plots
        gg.combined.n <- 
            n_plot_fun(
                combined_filter,
                "Combined filters",
                "black",
                x
            )
        
        gg.chimera.n <- 
            n_plot_fun(
                chimera_filter,
                "Chimera filter",
                "#9e0142",
                x
            )
        
        gg.length.n <- 
            n_plot_fun(
                length_filter,
                "Length filter",
                "#f46d43",
                x
            )
            
        gg.phmm.n <-
            n_plot_fun(
                phmm_filter,
                "PHMM filter",
                "#fee08b",
                x
            )
            
        gg.frame.n <- 
            n_plot_fun(
                frame_filter,
                "Frame/stop filter",
                "#abdda4",
                x
            )  
        
        # combine count plots
        
        gg.n_asvs <-
            gg.combined.n + gg.chimera.n + gg.length.n + gg.phmm.n + gg.frame.n +
            patchwork::plot_layout(ncol = 1, guides = "collect") +
            patchwork::plot_annotation(title = "Unique ASVs", subtitle = paste0("Flowcell: ", unique(x$read_group), "\nPCR primers: ", unique(x$primers)))
        
        ggsave(paste0("asv_count_",unique(x$read_group), "_", unique(x$primers),".pdf"), gg.n_asvs, width = 8, height = 12)

        }
    )

gc()

}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})