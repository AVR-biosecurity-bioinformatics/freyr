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

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "rt_samples",    
    "rt_group",
    "samplesheet_split_file" 
)
lapply(nf_vars, nf_var_check)

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

# get samplesheet_split as tibble
samplesheet_split <- readr::read_csv(samplesheet_split_file)

# add column names to tibble
colnames(sample_tibble) <- c("stage","sample_primers","read_group","primers","fwd","rev") 

sample_tibble_up <- 
    sample_tibble %>%
    # add "sample" from samplesheet
    dplyr::left_join(
      ., 
      samplesheet_split %>% dplyr::select(sample_primers, sample),
      by = "sample_primers"
    ) %>%
    dplyr::select(stage, sample, sample_primers, read_group, primers, fwd, rev) %>% 
    dplyr::arrange(sample_primers, primers, desc(fwd))

# reduce "input" (pre-split_loci) stage rows to one, changing pcr_primers to "combined" as reads have not been assigned by primer seq
sample_tibble_noinput <- 
    sample_tibble_up %>% 
    dplyr::filter(stage != "input") # remove "input" rows

sample_tibble_input <- 
    sample_tibble_up %>% 
    dplyr::filter(stage == "input") %>% 
    dplyr::group_by(sample) %>% 
    dplyr::slice(1) %>% # keep only one "input" row
    dplyr::ungroup() %>% 
    dplyr::mutate(
        sample_primers = sample, # change sample_primers to sample
        primers = "combined" # change pcr_primers as not used yet
    )

sample_tibble_altered <- rbind(sample_tibble_input, sample_tibble_noinput)

readr::write_csv(sample_tibble_altered, "sample_tibble.csv") # for debugging

## join group-level read tracking files into a single tibble

group_tibble <- tibble::tibble() # new tibble
for (i in 1:length(rt_group)) { # loop through .csv and add values to tibble as new rows
    new_csv <- readr::read_csv(rt_group[i], show_col_types = F)
    group_tibble <- rbind(group_tibble, new_csv)
}

group_tibble_up <- 
    group_tibble %>%
    # add "sample" from samplesheet
    dplyr::left_join(
      ., 
      samplesheet_split %>% dplyr::select(sample_primers, sample),
      by = "sample_primers"
    ) %>%
    dplyr::select(stage, sample, sample_primers, read_group, primers, pairs) %>% 
    dplyr::arrange(sample_primers, primers, desc(pairs))

readr::write_csv(group_tibble_up, "group_tibble.csv") # for debugging

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
    "filter_combined",
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

read_tracker_long <- 
    sample_tibble_altered %>%
    dplyr::select(-rev, pairs = fwd) %>% ### TODO: check fwd and rev counts are the same; for now assume and just use fwd as read pair count
    rbind(., group_tibble_up) 
    
read_tracker_inputcol <- 
    read_tracker_long %>% # pull out 
    dplyr::filter(primers == "combined") %>% 
    dplyr::select(sample, input = pairs)

read_tracker_wide <- 
    read_tracker_long %>% 
    dplyr::filter(primers != "combined") %>% 
    dplyr::left_join(., read_tracker_inputcol, by = "sample") %>% 
    tidyr::pivot_wider(names_from = stage, values_from = pairs) %>%
    dplyr::select(
      tidyselect::any_of(
        c(
          "sample",
          "sample_primers",
          "primers", 
          "read_group", 
          stage_vec
        )
      )
    )

readr::write_csv(read_tracker_wide, "read_tracker.csv")
readr::write_csv(read_tracker_long, "read_tracker_long.csv")

# create colour vector for primers
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
gg.read_tracker <- 
    read_tracker_long %>% 
    dplyr::mutate( 
        stage = forcats::fct_relevel(stage, stage_vec), # reorder x-axis
        primers = forcats::fct_relevel(primers, "combined") # make sure "combined" is always first
        ) %>% 
    ggplot2::ggplot(aes(x = stage, y = pairs, fill = primers)) +
    geom_col() + # TODO: Add % retention labels to the top of each bar
    scale_y_continuous(label = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1))) +
    scale_fill_manual(values = wong_pal) +
    facet_grid(read_group~.) +
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
        fill = "Primers"
    )

pdf(file="read_tracker.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker)
try(dev.off(), silent=TRUE)

### TODO: add a per-sample line graph output
gg.read_tracker_sample <- 
    read_tracker_wide %>% 
    tidyr::pivot_longer(
        cols = input:filter_sample_taxon,
        names_to = "stage",
        values_to = "pairs"
    ) %>% 
    dplyr::mutate( 
            stage = forcats::fct_relevel(stage, stage_vec), # reorder x-axis
            alpha = case_when( ### TODO: Use the locus parameter min reads as a filter here
                pairs < 1000 ~ 0.6,
                pairs >= 1000 ~ 1
                )
            ) %>% 
    ggplot2::ggplot(aes(x = stage, y = pairs, colour = primers)) +
    geom_line(aes(group = sample, alpha = alpha)) +
    scale_colour_manual(values = wong_pal[-1]) +
    scale_alpha_identity() +
    facet_grid(read_group~primers) +
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
        fill = "Primers"
        )

pdf(file="read_tracker_sample.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker_sample)
try(dev.off(), silent=TRUE)



### TODO: Make variant of per-sample line graph that is normalised by % of input reads per sample 
## and add average across flow cell and locus


# stop(" *** stopped manually *** ") ##########################################
