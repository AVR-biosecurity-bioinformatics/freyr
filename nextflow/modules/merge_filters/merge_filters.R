#!/usr/bin/env Rscript
.libPaths(.Library)

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

# Import functions
functions_file <- Sys.which("functions.R")
if (functions_file == "") {
    stop("Could not find functions.R on PATH") 
}
sys.source(functions_file, envir = .GlobalEnv)

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

### run R code

# combine seqtabs together (wide format)
seqtab_combined <- merge_count_matrices(
    lapply(
        stringr::str_split_1(seqtab_tibble_list, pattern = " "),
        function(f) {
            dt <- data.table::fread(f)
            seq_names <- dt$seq_name
            dt[, seq_name := NULL]
            mat <- as.matrix(dt)
            rownames(mat) <- seq_names
            mat
        }
    )
)

# Merge fasta together
merged_fasta <- merge_refseq(
    lapply(
        stringr::str_split_1(fasta_list, pattern = " "),
        Biostrings::readDNAStringSet
    )
)

# export combined .fasta of all sequences across flowcells 
write_fasta(
    merged_fasta,
    paste0(primers, "_combseqs.fasta")
)

# Create lookup table of sequence and names
seq_lookup <- data.table::data.table(
    seq_name = names(merged_fasta),
    sequence = as.character(merged_fasta)
)

# add sequence to seqtab and export
seqtab_export <- data.table::as.data.table(
    seqtab_combined,
    keep.rownames = "seq_name"
)

seqtab_export <- seq_lookup[
    seqtab_export,
    on = "seq_name"
]

data.table::fwrite(seqtab_export, paste0(primers, "_seqtab_combined.csv") )

# Combine filter outputs into one table with one row per sequence.
filters_combined <- data.table::rbindlist(
  lapply(
    stringr::str_split_1(filter_tibble_list, pattern = " "),
    data.table::fread
  ),
  fill = TRUE,
  use.names = TRUE
)[
  ,
  lapply(
    .SD,
    function(x) {
      x <- x[!is.na(x)]
      if (length(x)) x[1L] else NA
    }
  ),
  by = .(seq_name, sequence)
]

# Export the combined filter table.
data.table::fwrite(
  filters_combined,
  paste0(primers, "_filters.csv")
)

# Read sample metadata for the current primer set.
sample_read_group <- readr::read_csv(
  samplesheet_split,
  col_select = c(
    sample_primers,
    read_group,
    primers
  ),
  show_col_types = FALSE
) %>%
  dplyr::filter(primers == .env$primers) %>%
  dplyr::select(
    sample_primers,
    read_group
  ) %>%
  data.table::as.data.table()

# Combine sequence and filter metadata, retaining one row per ASV.
filter_metadata <- merge(
  seq_lookup,
  filters_combined,
  by = c("seq_name", "sequence"),
  all.x = TRUE,
  sort = FALSE
)

# Add the combined filter and static plotting metadata.
filter_metadata[
  ,
  `:=`(
    combined_filter =
      chimera_filter &
      length_filter &
      phmm_filter &
      frame_filter,
    length = nchar(sequence),
    primers = primers,
    concat = grepl("N{10}", sequence)
  )
]

# Align the metadata rows with the ASV rows in the abundance matrix.
filter_metadata <- filter_metadata[
  match(
    rownames(seqtab_combined),
    seq_name
  )
]

# Ensure every ASV in the abundance matrix has matching filter metadata.
if (anyNA(filter_metadata$seq_name)) {
  stop(
    "Some ASVs in seqtab_combined could not be matched ",
    "to filter_metadata."
  )
}

if (!identical(
  as.character(filter_metadata$seq_name),
  rownames(seqtab_combined)
)) {
  stop("filter_metadata is not aligned with seqtab_combined.")
}

# Ensure every sample in the metadata is present in the abundance matrix.
missing_samples <- setdiff(
  sample_read_group$sample_primers,
  colnames(seqtab_combined)
)

if (length(missing_samples)) {
  stop(
    "Samples are missing from seqtab_combined: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Used for lookup with .SD
filter_cols <- c(
  "chimera_filter",
  "length_filter",
  "phmm_filter",
  "frame_filter",
  "combined_filter"
)

# Convert logical filter results to ordered pass/fail factors once.
filter_metadata[
  ,
  (filter_cols) := lapply(
    .SD,
    function(x) {
      factor(
        ifelse(x, "pass", "fail"),
        levels = c("pass", "fail")
      )
    }
  ),
  .SDcols = filter_cols
]

read_groups <- unique(sample_read_group$read_group)

# Create one ASV summary table per read group.
filter_info <- lapply(
  read_groups,
  function(current_read_group) {

    # Get samples belonging to the current read group.
    current_samples <- sample_read_group[
      read_group == current_read_group,
      sample_primers
    ]

    # Calculate total abundance per ASV across the current samples.
    abundance <- rowSums(seqtab_combined[ , current_samples, drop = FALSE ])

    # Copy the shared metadata and add read-group-specific information.
    out <- data.table::copy(filter_metadata)
    out[, abundance := abundance]
    out[, read_group := current_read_group]

    # Keep observed ASVs and sort by decreasing abundance.
    out <- out[abundance > 0]
    data.table::setorder(out, -abundance)

    out
  }
)

# Name each list element by its read group.
names(filter_info) <- read_groups

# save filter summary files, one per flowcell
for (x in filter_info) {
    data.table::fwrite(x,paste0("asv_cleanup_", x$read_group[1], "_", x$primers[1],".csv"))
}

## export read tracking group tibble per flowcell (how many reads pass through each filter for each sample?)
filter_stages <- c(
  chimera_filter = "filter_chimera",
  length_filter = "filter_length",
  phmm_filter = "filter_phmm",
  frame_filter = "filter_frame",
  combined_filter = "filter_combined"
)

for (current_read_group in read_groups) {
  current_samples <- sample_read_group[ read_group == current_read_group, sample_primers]
  abundance <- seqtab_combined[, current_samples, drop = FALSE ]

  readsout <- data.table::rbindlist(
    lapply(
      names(filter_stages),
      function(filter_col) {

        passed <- filter_metadata[[filter_col]]

        data.table::data.table(
          stage = filter_stages[[filter_col]],
          sample_primers = current_samples,
          read_group = current_read_group,
          primers = primers,
          pairs = colSums(abundance[passed == "pass",,drop = FALSE])
        )
      }
    )
  )

  data.table::fwrite(readsout, paste0("filter_merged_", current_read_group, "_", primers, "_readsout.csv" )
  )
}

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