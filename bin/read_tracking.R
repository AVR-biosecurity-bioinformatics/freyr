#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "forcats",
    "ggplot2",
    "Hmisc",
    "readr",
    "scales",
    "stringr",
    "tibble",
    "tidyr",
    NULL
    )

invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

## check and define variables 
rt_samples <- # convert Groovy to R list format
    stringr::str_extract_all(rt_samples, pattern = "[^\\s,\\[\\]]+") %>% unlist()

rt_group <- # convert Groovy to R list format
    stringr::str_extract_all(rt_group, pattern = "[^\\s,\\[\\]]+") %>% unlist()

### run R code

## join sample-level read tracking files into a single tibble
sample_tibble <- tibble::tibble() # new tibble
for (i in 1:length(rt_samples)) { # loop through .csv and add values to tibble as new rows
    new_csv <- readr::read_csv(rt_samples[i], show_col_types = F, col_names = F)
    sample_tibble <- rbind(sample_tibble, new_csv)
}

colnames(sample_tibble) <- c("stage","sample_id","fcid","pcr_primers","fwd","rev") 

sample_tibble <- sample_tibble %>%
    dplyr::mutate(
        sample_id_com = sample_id, # sample_id before locus split, if done ("com" for "combined")
        sample_id = paste0(sample_id,"_",pcr_primers) # make sample_id consistent with "sample_id" + "pcr_primers" format
        ) %>% 
    dplyr::select(stage, sample_id_com, sample_id, fcid, pcr_primers, fwd, rev) %>% 
    dplyr::arrange(sample_id, pcr_primers, desc(fwd))

sample_id_matching <- sample_tibble %>% dplyr::select(sample_id_com, sample_id) %>% dplyr::distinct() # get tibble of sample_id and matching sample_id_com

# reduce "input" (pre-split_loci) stage rows to one, changing pcr_primers to "combined" as reads have not been assigned by primer seq
sample_tibble_noinput <- sample_tibble %>% 
    dplyr::filter(stage != "input") # remove "input" rows

sample_tibble_input <- sample_tibble %>% 
    dplyr::filter(stage == "input") %>% 
    dplyr::group_by(sample_id_com) %>% 
    dplyr::slice(1) %>% # keep only one "input" row
    dplyr::ungroup() %>% 
    dplyr::mutate(
        sample_id = sample_id_com, # change sample_id as pcr_primers not used
        pcr_primers = "combined" # change pcr_primers as not used yet
    )

sample_tibble_altered <- rbind(sample_tibble_input, sample_tibble_noinput)

readr::write_csv(sample_tibble_altered, "sample_tibble.csv") # for debugging

## join group-level read tracking files into a single tibble

group_tibble <- tibble::tibble() # new tibble
for (i in 1:length(rt_group)) { # loop through .csv and add values to tibble as new rows
    new_csv <- readr::read_csv(rt_group[i], show_col_types = F)
    group_tibble <- rbind(group_tibble, new_csv)
}

group_tibble <- group_tibble %>%
    dplyr::mutate(
        sample_id = ifelse(
            stringr::str_detect(sample_id, pcr_primers), # if sample_id contains pcr_primers...
            sample_id, # keep same
            paste0(sample_id,"_",pcr_primers) # else add pcr_primers to the end of sample_id
            )
        ) %>% # make sample_id consistent with "sample_id" + "pcr_primers" format
    dplyr::left_join(., sample_id_matching, by = "sample_id") %>% # add sample_id_com to tibble
    dplyr::select(stage, sample_id_com, sample_id, fcid, pcr_primers, pairs) %>% 
    dplyr::arrange(sample_id, pcr_primers, desc(pairs))

readr::write_csv(group_tibble, "group_tibble.csv") # for debugging

## combine sample and group tibbles together
stage_vec <- c(
    "input", 
    "split_loci", 
    "primer_trim",
    "read_filter",
    "dada_mergereads", 
    "filter_chimera", 
    "filter_length", 
    "filter_phmm", 
    "filter_frame", 
    "filter_seqtab",
    "classified_root", 
    "classified_kingdom", 
    "classified_phylum",
    "classified_class",
    "classified_order", 
    "classified_family", 
    "classified_genus", 
    "classified_species", 
    "filter_sample_taxon"
    )

read_tracker_long <- sample_tibble_altered %>%
    dplyr::select(-rev, pairs = fwd) %>% ### TODO: check fwd and rev counts are the same; for now assume and just use fwd as read pair count
    rbind(., group_tibble) 
    
read_tracker_inputcol <- read_tracker_long %>% # pull out 
    dplyr::filter(pcr_primers == "combined") %>% 
    dplyr::select(sample_id_com, input = pairs)

read_tracker_wide <- read_tracker_long %>% 
    dplyr::filter(pcr_primers != "combined") %>% 
    dplyr::left_join(., read_tracker_inputcol, by = "sample_id_com") %>% 
    tidyr::pivot_wider(names_from = stage, values_from = pairs) %>%
    dplyr::select(tidyselect::any_of(c(
        "sample_id_com",
        "sample_id",
        "pcr_primers", 
        "fcid", 
        stage_vec
    )))

readr::write_csv(read_tracker_wide, "read_tracker.csv")
readr::write_csv(read_tracker_long, "read_tracker_long.csv")

# create colour vector for pcr_primers
wong_pal <- c(
    "grey40", # grey
    "#E69F00", # orange
    "#56B4E9", # light blue
    "#009E73", # green
    "#F0E442", # yellow
    "#0072B2", # dark blue
    "#D55E00", # red-orange
    "#CC79A7" # pink
    )

## plot read tracking by locus
gg.read_tracker <- read_tracker_long %>% 
    dplyr::mutate( 
        stage = forcats::fct_relevel(stage, stage_vec), # reorder x-axis
        pcr_primers = forcats::fct_relevel(pcr_primers, "combined") # make sure "combined" is always first
        ) %>% 
    ggplot2::ggplot(aes(x = stage, y = pairs, fill=pcr_primers)) +
    geom_col() + # TODO: Add % retention labels to the top of each bar
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(values = wong_pal) +
    facet_grid(fcid~.) +
    theme_bw() +
    theme(
        strip.background = element_rect(colour = "black", fill = "lightgray"),
        strip.text = element_text(size=9, family = ""),
        axis.text.x =element_text(angle=45, hjust=1, vjust=1),
        plot.background = element_blank(),
        text = element_text(size=9, family = ""),
        axis.text = element_text(size=8, family = ""),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid = element_line(linewidth = rel(0.5)),
    ) +
    labs(
        x = "Pipeline step",
        y = "Reads retained",
        fill = "PCR primers"
        )

pdf(file="read_tracker.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker)
try(dev.off(), silent=TRUE)

### TODO: add a per-sample line graph output
gg.read_tracker_sample <- read_tracker_wide %>% 
    tidyr::pivot_longer(
        cols = input:filter_sample_taxon,
        names_to = "stage",
        values_to = "pairs"
    ) %>% 
    dplyr::mutate( 
            stage = forcats::fct_relevel(stage, stage_vec), # reorder x-axis
            alpha = case_when( ### TODO: Use the locus parameter min reads as a filter here
                pairs < 1000 ~ 0.8,
                pairs >= 1000 ~ 1
                )
            ) %>% 
    ggplot2::ggplot(aes(x = stage, y = pairs, colour = pcr_primers)) +
    geom_line(aes(group = sample_id_com, alpha = alpha)) +
    scale_colour_manual(values = wong_pal[-1]) +
    facet_grid(fcid~pcr_primers) +
    theme_bw() +
    theme(
        strip.background = element_rect(colour = "black", fill = "lightgray"),
        strip.text = element_text(size=9, family = ""),
        axis.text.x =element_text(angle=45, hjust=1, vjust=1),
        plot.background = element_blank(),
        text = element_text(size=9, family = ""),
        axis.text = element_text(size=8, family = ""),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.grid = element_line(linewidth = rel(0.5)),
    ) +
    labs(
        x = "Pipeline step",
        y = "Reads retained",
        fill = "PCR primers"
        )

pdf(file="read_tracker_sample.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker_sample)
try(dev.off(), silent=TRUE)

### TODO: Make variant of per-sample line graph that is normalised by % of input reads per sample 
## and add average across flow cell and locus


# stop(" *** stopped manually *** ") ##########################################
