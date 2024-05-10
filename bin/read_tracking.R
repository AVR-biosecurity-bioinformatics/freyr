#!/usr/bin/env Rscript

# code produces:
# - plot of N reads retained after each step, faceted by flowcell, stacked bars coloured by locus
# - .csv file that show same info by per sample (split by locus)

### required inputs:
# - unfiltered phyloseq
# - filtered phyloseq
# - reads input for trimming
# - 


## check and define variables 
rt_samples <- # convert Groovy to R list format
    stringr::str_extract_all(rt_samples, pattern = "[^\\s,\\[\\]]+") %>% unlist()

rt_group <- # convert Groovy to R list format
    stringr::str_extract_all(rt_group, pattern = "[^\\s,\\[\\]]+") %>% unlist()

### run R code

## join sample-level read tracking files into a single tibble
sample_tibble <- tibble() # new tibble
for (i in 1:length(rt_samples)) { # loop through .csv and add values to tibble as new rows
    new_csv <- read_csv(rt_samples[i], show_col_types = F, col_names = F)
    sample_tibble <- rbind(sample_tibble, new_csv)
}

colnames(sample_tibble) <- c("stage","sample_id","fcid","pcr_primers","fwd","rev") 

sample_tibble <- sample_tibble %>%
    mutate(
        sample_id_com = sample_id, # sample_id before locus split, if done ("com" for "combined")
        sample_id = paste0(sample_id,"_",pcr_primers) # make sample_id consistent with "sample_id" + "pcr_primers" format
        ) %>% 
    dplyr::select(stage, sample_id_com, sample_id, fcid, pcr_primers, fwd, rev) %>% 
    dplyr::arrange(sample_id, pcr_primers, desc(fwd))

sample_id_matching <- sample_tibble %>% dplyr::select(sample_id_com, sample_id) %>% distinct() # get tibble of sample_id and matching sample_id_com

# reduce "input" (pre-split_loci) stage rows to one, changing pcr_primers to "combined" as reads have not been assigned by primer seq
sample_tibble_noinput <- sample_tibble %>% 
    dplyr::filter(stage != "input") # remove "input" rows

sample_tibble_input <- sample_tibble %>% 
    dplyr::filter(stage == "input") %>% 
    group_by(sample_id_com) %>% 
    dplyr::slice(1) %>% # keep only one "input" row
    ungroup() %>% 
    dplyr::mutate(
        sample_id = sample_id_com, # change sample_id as pcr_primers not used
        pcr_primers = "combined" # change pcr_primers as not used yet
    )

sample_tibble_altered <- rbind(sample_tibble_input, sample_tibble_noinput)

write_csv(sample_tibble_altered, "sample_tibble.csv") # for debugging

## join group-level read tracking files into a single tibble

group_tibble <- tibble() # new tibble
for (i in 1:length(rt_group)) { # loop through .csv and add values to tibble as new rows
    new_csv <- read_csv(rt_group[i], show_col_types = F)
    group_tibble <- rbind(group_tibble, new_csv)
}

group_tibble <- group_tibble %>%
    mutate(
        sample_id = ifelse(
            stringr::str_detect(sample_id, pcr_primers), # if sample_id contains pcr_primers...
            sample_id, # keep same
            paste0(sample_id,"_",pcr_primers) # else add pcr_primers to the end of sample_id
            )
        ) %>% # make sample_id consistent with "sample_id" + "pcr_primers" format
    left_join(., sample_id_matching, by = "sample_id") %>% # add sample_id_com to tibble
    dplyr::select(stage, sample_id_com, sample_id, fcid, pcr_primers, pairs) %>% 
    dplyr::arrange(sample_id, pcr_primers, desc(pairs))

write_csv(group_tibble, "group_tibble.csv") # for debugging

## combine sample and group tibbles together
steps_vec <- c(
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
    
read_tracker_inputcol <- read_tracker_long %>%
    dplyr::filter(pcr_primers == "combined") %>% 
    dplyr::select(sample_id_com, input_reads = pairs)

read_tracker_wide <- read_tracker_long %>% 
    dplyr::filter(pcr_primers != "combined") %>% 
    left_join(., read_tracker_inputcol, by = "sample_id_com") %>% 
    pivot_wider(names_from = stage, values_from = pairs) %>%
    dplyr::select(any_of(c(
        "sample_id_com",
        "sample_id",
        "pcr_primers", 
        "fcid", 
        "input_reads",
        steps_vec
    )))

write_csv(read_tracker_wide, "read_tracker.csv")

## plot read tracking
gg.read_tracker <- read_tracker_long %>% 
    dplyr::mutate(step = factor(step, levels=steps_vec)) %>% # reorder step factor
    ggplot(aes(x = step, y = pairs, fill=pcr_primers))+
    geom_col() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
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
    labs(x = "Pipeline step",
            y = "Reads retained",
            fill = "PCR primers")
    pdf(file="read_tracker.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker)
    try(dev.off(), silent=TRUE)

stop(" *** stopped manually *** ") ##########################################
