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
        sample_id = paste0(sample_id,"_",pcr_primers)) %>% # make sample_id consistent with "sample_id" + "pcr_primers" format
    dplyr::arrange(sample_id, pcr_primers, desc(fwd))

sample_id_matching <- sample_tibble %>% dplyr::select(sample_id_com, sample_id) %>% distinct() # get tibble of sample_id and matching sample_id_com

write_csv(sample_tibble, "sample_tibble.csv") # for debugging

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

read_tracker <- sample_tibble %>%
    dplyr::select(-rev, pairs = fwd) %>% ### TODO: check fwd and rev counts are the same; for now assume and just use fwd as read pair count
    rbind(., group_tibble) %>%
    pivot_wider(names_from = stage, values_from = pairs) %>%
    dplyr::select(any_of(c(
        "sample_id",
        "pcr_primers", 
        "fcid", 
        steps_vec
    )))

write_csv(read_tracker, "read_tracker.csv")

## plot read tracking
gg.read_tracker <- read_tracker %>%
    pivot_longer(cols = -c("sample_id", "fcid", "pcr_primers"),
                    names_to = "step",
                    values_to= "reads") %>%
    group_by(sample_id) %>%
    # group_modify(~{ # don't need this step because sample_id should be unique
    #     # When a sample name shares multiple sample ids, select a single sample to avoid double counting
    #     if(length(unique(.x$pcr_primers)) >1){
    #     .x %>%
    #         dplyr::filter(!step == "input_reads") %>%
    #         bind_rows(.x %>%
    #                     dplyr::filter(step == "input_reads") %>%
    #                     mutate(pcr_primers="Mixed",
    #                             sample_id=NA_character_) %>%
    #                     dplyr::slice(1))
    #     } else {
    #     .x
    #     }
    # }) %>%
    dplyr::mutate(step = factor(step, levels=steps_vec)) %>% # reorder step factor
    ggplot(aes(x = step, y = reads, fill=pcr_primers))+
    geom_col() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))+
    facet_grid(fcid~.)+
    theme_bw()+
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
