#### NOTE:
#### The following variables are defined in the module script blocks:
#### projectDir : the project directory absolute path; location of main.nf and "working directory"
#### data_loc : the name of the data directory (usually "data" but could be "test_data")


#### jack's functions

## checks existence of a variable and print error message if not defined
nf_var_check <- function(x) {
  if (!exists(x)) {
    stop(paste0("The variable'",x,"' is not defined! Make sure to check the Nextflow process inputs."))
  } else {
    print(paste0("Input variable '",x,"' = ",eval(parse(text = x))))
  }
}

## collapses repetive Groovy list variables down to a single variable
parse_nf_var_repeat <- function(x) {
  variable <- 
    stringr::str_extract_all(
      x, 
      pattern = "[^\\s,\\[\\]]+" # extract all runs of characters that aren't ' ' ',' '[' or ']' 
    ) %>% 
    unlist() %>%
    tibble::as_tibble_col(column_name = "col") %>% 
    unique() %>%
    dplyr::pull(col)
  
  if (length(variable) == 1) {
    out <- variable
  } else {
    out <- stop("*** nf variable contains multiple unique values! ***")
  }
  return(out)
}

# Coalesce taxonomy function
coalesce_tax <- function (x, y, suffix = c(".x", ".y"), prefer="left", join = dplyr::full_join, 
                          ...) {
  if(!"OTU" %in% colnames(x)){
    x <- x %>%
      tibble::as_tibble(rownames = "OTU")
  }
  if(!"OTU" %in% colnames(y)){
    y <- y %>%
      tibble::as_tibble(rownames = "OTU")
  }
  joined <- join(x, y, by = "OTU", suffix = suffix, ...)
  cols <- union(names(x), names(y))
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  to_coalesce <- unique(substr(to_coalesce, 1, nchar(to_coalesce) - nchar(suffix_used)))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~{
    left_side <- joined[[paste0(.x, suffix[1])]]
    right_side <- joined[[paste0(.x, suffix[2])]]
    if ((!all(is.na(left_side)) && !all(is.na(right_side))) && 
        (!class(left_side) == class(right_side))) {
      class(right_side) <- class(left_side)
    }
    # Set the right side to NA's if they are present in the left side
    if(prefer == "left"){
      right_side[!is.na(left_side)] <- NA
    }else if (prefer == "right"){
      left_side[!is.na(right_side)] <- NA
    }
    dplyr::coalesce(right_side, left_side)
  })
  names(coalesced) <- to_coalesce
  out <- dplyr::bind_cols(joined, coalesced)[cols] %>%
    tibble::column_to_rownames("OTU")
  return(out)
}

# Safely bind two data frames/tibbles that may contain the same columns with different classes.
new_bind <- function(a, b) {
  common_cols <- intersect(names(a), names(b))
  b[common_cols] <- map2_df(b[common_cols], map(a[common_cols], class), ~{class(.x) <- .y;.x})
  bind_rows(a, b)  
}

# New subset taxa function that allows variable column inputs
subset_taxa_new <- function(physeq, rank, value){
  if (is.null(phyloseq::tax_table(physeq))) {
    cat("Nothing subset. No taxonomyTable in physeq.\n")
    return(physeq)
  } else {
    oldMA <- as(phyloseq::tax_table(physeq), "matrix")
    oldDF <- data.frame(oldMA)
    newDF <- oldDF %>%
      dplyr::filter(UQ(sym(rank)) == value)
    newMA <- as(newDF, "matrix")
    phyloseq::tax_table(physeq) <- phyloseq::tax_table(newMA)
    return(physeq)
  }
}

# New merge phyloseq function that accepts a list of phyloseq objects
merge_phyloseq_new <- function(arguments) {
  comp.list <- unlist(
    lapply(arguments, phyloseq:::splat.phyloseq.objects),
    recursive = FALSE
  )
  comp.groups <- split(comp.list, names(comp.list))
  merged.list <- lapply(comp.groups, function(x) {
    
    if (length(x) == 1L) {
      x[[1]]
    } else {
      Reduce(
        phyloseq::merge_phyloseq_pair,
        x
      )
    }
  })
  do.call(phyloseq::phyloseq, unname(merged.list))
}

# Merge reference sequences, retaining the first sequence for each ASV name.
merge_refseq <- function(x) {
  
  seqs <- do.call(
    c,
    lapply(x, function(z) {
      out <- as.character(z)
      names(out) <- names(z)
      out
    })
  )
  seqs <- seqs[!duplicated(names(seqs))]
  names(seqs) <- sub("^refseq\\.", "", names(seqs)) # remove refseq prefix
  Biostrings::DNAStringSet(seqs)
}

# Merge integer count matrices with taxa as rows and disjoint samples.
merge_count_matrices <- function(x) {
  if (!length(x)) {
    stop("No count matrices were provided.")
  }

  taxa <- unique(unlist(
    lapply(x, rownames),
    use.names = FALSE
  ))

  samples <- unlist(
    lapply(x, colnames),
    use.names = FALSE
  )

  taxa_index <- stats::setNames(
    seq_along(taxa),
    taxa
  )

  # Allocate the final abundance matrix once.
  merged <- matrix(
    0L,
    nrow = length(taxa),
    ncol = length(samples),
    dimnames = list(taxa, samples)
  )

  column_start <- 1L

  for (i in seq_along(x)) {
    current <- x[[i]]

    if (is.null(rownames(current))) {
      stop("Count matrix ", i, " has no taxon row names.")
    }

    if (is.null(colnames(current))) {
      stop("Count matrix ", i, " has no sample column names.")
    }

    # Convert legacy double-backed count matrices if necessary.
    if (typeof(current) != "integer") {
      storage.mode(current) <- "integer"
    }

    row_index <- unname(
      taxa_index[rownames(current)]
    )

    if (anyNA(row_index)) {
      stop(
        "Failed to match all taxa from count matrix ",
        i,
        " to the merged taxon set."
      )
    }

    column_end <- column_start + ncol(current) - 1L

    merged[
      row_index,
      seq.int(column_start, column_end)
    ] <- current

    column_start <- column_end + 1L
  }

  merged
}


# Merge phyloseq function for disjoint samples that accepts a list of phyloseq objects
merge_phyloseq_disjoint_samples <- function(arguments) {
  if (!length(arguments)) {
    stop("arguments must contain at least one object.")
  }
  
  arguments <- unname(arguments)
  
  # Require sample names to be unique across all input objects.
  samples <- unlist(
    lapply(arguments, phyloseq::sample_names),
    use.names = FALSE
  )
  
  duplicated_samples <- unique(
    samples[duplicated(samples)]
  )
  
  if (length(duplicated_samples)) {
    stop(
      "Sample names must not overlap between input phyloseq objects. ",
      "Duplicated sample names: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }
  
  # Merge integer OTU tables with taxa as rows and disjoint samples.
  merge_otu_tables <- function(x) {
    taxa_rows <- vapply(
      x,
      phyloseq::taxa_are_rows,
      logical(1)
    )

    if (!all(taxa_rows)) {
      stop(
        "All OTU tables must have taxa as rows. ",
        "Incorrect input indices: ",
        paste(which(!taxa_rows), collapse = ", ")
      )
    }

    matrices <- lapply(
      x,
      function(current) {
        as(current, "matrix")
      }
    )

    phyloseq::otu_table(
      merge_count_matrices(matrices),
      taxa_are_rows = TRUE
    )
  }
  
  # Merge sample metadata with disjoint sample rows.
  merge_sample_data <- function(x) {
    merged <- data.table::rbindlist(
      lapply(
        x,
        function(current) {
          data.table::as.data.table(
            as(current, "data.frame"),
            keep.rownames = "sample_name"
          )
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
    
    sample_ids <- merged[["sample_name"]]
    merged[, sample_name := NULL]
    
    merged <- as.data.frame(merged)
    rownames(merged) <- sample_ids
    
    phyloseq::sample_data(merged)
  }
  
  # Merge taxonomy tables, retaining the first assignment for duplicate ASVs.
  # Assumes taxonomy assignments for the same ASV are consistent across inputs.
  merge_tax_tables <- function(x) {
    merged <- data.table::rbindlist(
      lapply(
        x,
        function(current) {
          data.table::as.data.table(
            as(current, "matrix"),
            keep.rownames = "seq_name"
          )
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
    
    merged <- merged[!duplicated(seq_name)]
    
    seq_names <- merged[["seq_name"]]
    merged[, seq_name := NULL]
    
    merged <- as.matrix(merged)
    rownames(merged) <- seq_names
    
    phyloseq::tax_table(merged)
  }
    
  # Split each phyloseq object into its constituent components.
  comp.list <- unlist(
    lapply(
      arguments,
      phyloseq:::splat.phyloseq.objects
    ),
    recursive = FALSE
  )
  
  comp.groups <- split(
    comp.list,
    names(comp.list)
  )
  
  merged.list <- lapply(
    names(comp.groups),
    function(component) {
      x <- comp.groups[[component]]
      
      if (length(x) == 1L) {
        return(x[[1L]])
      }
      
      if (component == "otu_table") {
        merge_otu_tables(x)
      } else if (component == "sam_data") {
        merge_sample_data(x)
      } else if (component == "tax_table") {
        merge_tax_tables(x)
      } else if (component == "refseq") {
        merge_refseq(x)
      } else {
        # Retain standard phyloseq behaviour for other components,
        # such as phylogenetic trees.
        Reduce(
          phyloseq::merge_phyloseq_pair,
          x
        )
      }
    }
  )
  
  do.call(
    phyloseq::phyloseq,
    unname(merged.list)
  )
}

summarise_phyloseq <- function(ps) {
  otutab <- as(phyloseq::otu_table(ps), "matrix")
  storage.mode(otutab) <- "integer" # Force integer storage of discrete count matrix
  if (!phyloseq::taxa_are_rows(ps)) otutab <- t(otutab)
  
  sample_cols <- colnames(otutab)
  
  df <- data.table::as.data.table(
    otutab,
    keep.rownames = "seq_name"
  )
  
  df[, sequence := as.character(phyloseq::refseq(ps)[seq_name])]
  
  taxtab <- data.table::as.data.table(
    as(phyloseq::tax_table(ps), "matrix"),
    keep.rownames = "seq_name"
  )
  
  tax_cols <- setdiff(names(taxtab), "seq_name")
  
  df <- taxtab[df, on = "seq_name"]
  
  data.table::setcolorder(
    df,
    c("seq_name", "sequence", tax_cols, sample_cols)
  )
  
  df
}

# Convert a phyloseq object into a long data.table containing one row per non-zero sample–taxon observation.
melt_phyloseq <- function(ps) {
  
  ## Abundance table-
  otutab <- as(phyloseq::otu_table(ps), "matrix")
  storage.mode(otutab) <- "integer" # Force integer storage of discrete count matrix
  if(!phyloseq::taxa_are_rows(ps)){otutab <- t(otutab)}
  
  # Identify only non-zero observations.
  index <- which(otutab > 0, arr.ind = TRUE)
  
  # If no non-zero obs, exit early
  if (length(index) == 0L) {return(data.table::data.table()) }
  
  # Convert otu table to tall form (one sample-taxon per row)
  df <- data.table::data.table(
    seq_name = rownames(otutab)[index[, "row"]],
    sample_primers = colnames(otutab)[index[, "col"]],
    Abundance = otutab[index]
  )
  
  ## Add Sample metadata
  if (length(phyloseq::sample_variables(ps)) > 0L) {
    samdf <- data.table::as.data.table(
      as(phyloseq::sample_data(ps), "data.frame"),
      keep.rownames = ".sample_name"
    )
    
    # Match metadata against sample names in the abundance table.
    sample_index <- match(df$sample_primers, samdf$.sample_name)
    
    # Avoid duplicating the sample identifier if sample_data already has sample_primers
    sample_columns <- setdiff(names(samdf), c(".sample_name", "sample_primers"))
    if (length(sample_columns) > 0L) {
      df <- cbind(
        df,
        samdf[sample_index, sample_columns, with = FALSE]
      )
    }
  }
  
  ## Add Taxonomy table
  if (length(phyloseq::rank_names(ps)) > 0L) {
    taxtab <- data.table::as.data.table(
      as(phyloseq::tax_table(ps), "matrix"),
      keep.rownames = ".seq_name"
    )
    taxon_index <- match(df$seq_name, taxtab$.seq_name)
    taxonomy_columns <- setdiff(names(taxtab),".seq_name")
    
    if (length(taxonomy_columns) > 0L) {
      df <- cbind(df, taxtab[taxon_index, taxonomy_columns, with = FALSE] )
    }
  }
  
  ## Add sequences
  seqs <- phyloseq::refseq(ps)
  if (!is.null(seqs)) {
    sequence_index <- match(df$seq_name, names(seqs))
    df[, sequence := as.character(seqs[sequence_index])]
  }
  
  # Arrange by Abundance, then sequence names (to approx. phyloseq behavior)
  data.table::setorder(df, -Abundance, seq_name )
  return(df)
}

# Write raw phyloseq
write_raw_phyloseq <- function(
    ps,
    cluster_dt,
    filename
) {
  
  raw_dt <- cluster_dt[
    melt_phyloseq(ps),
    on = "seq_name"
  ]
  
  data.table::setcolorder(
    raw_dt,
    c(
      "seq_name",
      "sequence",
      "cluster",
      setdiff(
        names(raw_dt),
        c(
          "seq_name",
          "sequence",
          "cluster"
        )
      )
    )
  )
  
  data.table::fwrite(
    raw_dt,
    filename
  )
  
  invisible(NULL)
}

# Write phyloseq summary
write_phyloseq_summary <- function(
    ps,
    cluster_dt,
    filename,
    tax_cols,
    filters = NULL
) {
  summary_dt <- summarise_phyloseq(ps)
  
  # Optionally add sequence-filter results.
  if (!is.null(filters)) {
    summary_dt <- filters[
      summary_dt,
      on = c("seq_name", "sequence")
    ]
  }
  
  # Add cluster assignments.
  summary_dt <- cluster_dt[
    summary_dt,
    on = "seq_name"
  ]
  
  preferred_cols <- c(
    "seq_name",
    "sequence",
    "cluster",
    tax_cols,
    "chimera_filter",
    "length_filter",
    "phmm_filter",
    "frame_filter"
  )
  
  preferred_cols <- intersect(
    preferred_cols,
    names(summary_dt)
  )
  
  data.table::setcolorder(
    summary_dt,
    c(
      preferred_cols,
      setdiff(names(summary_dt), preferred_cols)
    )
  )
  
  data.table::fwrite(
    summary_dt,
    filename
  )
  
  invisible(NULL)
}

# Rarefaction / species accumulation plot
rareplot <- function(ps, step="auto", threshold=0){
  if(step == "auto"){
    step <- ceiling(median(phyloseq::sample_sums(ps)) / 100)
  } else if (is.numeric(step)){
    step <- as.numeric(ceiling(step))
  } else {
    stop("Step must be an integer or 'auto' ")
  }
  if(phyloseq::taxa_are_rows(ps)){
    otutab <- t(as(phyloseq::otu_table(ps), "matrix"))
    storage.mode(otutab) <- "integer" # Force integer storage of discrete count matrix
  } else {
    otutab <- as(phyloseq::otu_table(ps), "matrix")
    storage.mode(otutab) <- "integer" # Force integer storage of discrete count matrix
  }
  otutab <- otutab[rowSums(otutab) > 0, colSums(otutab) > 0]
  rare <- vegan::rarecurve(otutab, step=step, tidy=TRUE)%>%
    dplyr::rename(sample_primers = Site, count = Sample, ASV = Species) %>%
    dplyr::left_join(
      data.frame(sample_primers = phyloseq::sample_data(ps)$sample_primers,
                 read_group = phyloseq::sample_data(ps)$read_group),
      by = "sample_primers"
    )
  
  gg.rare <- rare %>%
    ggplot2::ggplot() +
    geom_line(aes(x = count, y = ASV, group=sample_primers), alpha=0.3)+
    geom_point(data = rare %>% 
                 group_by(sample_primers) %>% 
                 top_n(1, count),
               aes(x = count, y = ASV, colour=count > threshold)) +
    scale_x_continuous(label = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))+
    scale_colour_manual(values=c("FALSE" = "#F8766D", "TRUE"="#619CFF"))+
    facet_wrap(read_group~., scales="free", ncol=1)+
    theme_bw()+
    theme(
      strip.background = element_rect(colour = "black", fill = "lightgray"),
      strip.text = element_text(size=9, family = ""),
      plot.background = element_blank(),
      text = element_text(size=9, family = ""),
      axis.text = element_text(size=8, family = ""),
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.grid = element_line(size = rel(0.5)),
    ) + labs(x = "Sequence reads",
             y = "Observed ASVs",
             colour = "Above sample filtering theshold") 
  
  return(gg.rare)
}

#' write DNAbin or AAbin to fasta
write_fasta <- function(x, file = "", compress = FALSE) {
  if(stringr::str_detect(file, "\\.gz$") & !compress){
    compress <- TRUE
  }
  if (!is.null(dim(x))) {
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if (inherits(x, "DNAbin")) {
    tmp <- DNAbin2char(x)
  } else if (is.list(x)) {
    if (length(x[[1]]) == 1L) {
      tmp <- unlist(x, use.names = TRUE)
    }
    else {
      tmp <- sapply(x, paste0, collapse = "")
    }
  } else {
    tmp <- x
  }
  reslen <- 2 * length(tmp)
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0(">", names(tmp))
  res[seq(2, reslen, by = 2)] <- tmp
  if(!file == ""){
    f <- if(compress){
      gzfile(file, "w")
    } else {
      file(file, "w")
    }
    writeLines(res, f)
    close(f)
  } else {
    writeLines(res)
  }
  invisible(NULL)
}

# Summarise_index ---------------------------------------------------------

#' Summarise index
#'
#' @param fq A path to a fastq file
#' @param qualityType The quality encoding of the fastq file
#'
#' @return
#' @export
#' @import ShortRead
#' @import tidyr
#' @import dplyr
#' @import stringr
#'
#' @examples
summarise_index <- function(fq, qualityType = "Auto") {
  fF <- ShortRead::FastqStreamer(fq)
  on.exit(close(fF))
  
  #Stream through fasta
  while( length(suppressWarnings(fqF <- yield(fF, qualityType = qualityType)))) {
    idF <- ShortRead::id(fqF)
  }
  
  # Check if index is present in fastq header
  index_check <- stringr::str_extract(as.character(idF[[1]]), pattern="(?!:)(?:.(?!:))+$")
  if(stringr::str_detect(index_check, "[A-Z]+[A-Z]")){
    index <- stringr::str_extract(as.character(idF), pattern="(?!:)(?:.(?!:))+$") %>%
      table() %>%
      as.data.frame() %>%
      tidyr::separate(col=1, into=c("index", "index2")) %>%
      dplyr::arrange(desc(Freq))
  } else{
    warning(paste0("Index sequences were not present in fastq header for sample: ", fq))
    index <- data.frame(index = NA_character_, index2=NA_character_, Freq = NA_integer_)
  }
  return(index)
}

# Get Qual Stats ----------------------------------------------------------

#' Get qual stats
#'
#' @param input
#' @param n The number of reads to sample
#'
#' @return
#' @export
#'
#' @examples
get_qual_stats <- function (input, n = 5e+05) {
  
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0),
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0),
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0),
                      rclabel = character(0), rc = numeric(0), file = character(0))
  #Check files arent empty
  fl <- input[file.size(input) > 28]
  if (length(input) > length(fl)) {
    message("Warning, the following files were empty: \n")
    message(paste0(input[!input %in% fl], " \n"))
  }
  # Start loop
  FIRST <- TRUE
  for (f in fl[!is.na(fl)]) {
    fq <- FastqFile(f)
    rc <- countFastq(f)$records
    close(fq)
    coll <- QACollate(QAFastqSource(f, n=n), QAQualityByCycle())
    srqa <- qa2(coll, BPPARAM=SerialParam(), verbose=FALSE)
    df <- data.frame(srqa$QAQualityByCycle@values)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }  else {      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count,  df$Cycle)
    
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(x) get_quant(x$Score, x$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(x) sum(x$Count), simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    } else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }
    
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)),
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s),
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc,
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score),
                                     label = basename(f), rclabel = rclabel, rc = rc,
                                     file = basename(f)))
  }
  
  anndf$minScore <- min(anndf$minScore)
  
  get_ee <- function(Q){
    ee <- 10^(-Q/ 10)
    return(ee)
  }
  
  
  plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
  plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
  means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
  EEmeans <- get_ee(rowsum(plotdf.summary$Score * plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle))
  q10s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.1), simplify = TRUE)
  EE10s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.1)), simplify = TRUE)
  q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.25), simplify = TRUE)
  EE25s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.25)), simplify = TRUE)
  q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.5), simplify = TRUE)
  EE50s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.50)), simplify = TRUE)
  q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.75), simplify = TRUE)
  EE75s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.75)), simplify = TRUE)
  q90s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_quant(x$Score, x$Count, 0.9), simplify = TRUE)
  EE90s <- by(plotdf.summary, plotdf.summary$Cycle, function(x) get_ee(get_quant(x$Score, x$Count, 0.9)), simplify = TRUE)
  cums <- by(plotdf.summary, plotdf.summary$Cycle, function(x) sum(x$Count), simplify = TRUE)
  statdf.summary <- data.frame(Cycle = as.integer(rownames(means)),
                               QMean = means,
                               EEmean = EEmeans,
                               Q10 = as.vector(q10s),
                               EE10 = as.vector(EE10s),
                               Q25 = as.vector(q25s),
                               EE25 = as.vector(EE25s),
                               Q50 = as.vector(q50s),
                               EE50 = as.vector(EE50s),
                               Q75 = as.vector(q75s),
                               EE75 = as.vector(EE75s),
                               Q90 = as.vector(q90s),
                               EE90 = as.vector(EE90s),
                               reads = 10 * as.vector(cums)/rc)
  attr(statdf.summary, "qsummary") <- "qsummary"
  return(statdf.summary)
}



# train IDTAXA -----------------------------------------------------------

train_idtaxa <- function(
    x,
    max_group_size = 10,
    max_iterations = 3,
    allow_group_removal = TRUE,
    orient = FALSE,
    drop_malformed = FALSE,
    quiet = FALSE
) {
  time <- Sys.time()
  
  # Convert input
  if (methods::is(x, "DNAbin")) {
    seqs <- DNAbin2DNAstringset(x)
  } else if (methods::is(x, "DNAStringSet")) {
    seqs <- x
  } else {
    stop("'x' must be a DNAbin or DNAStringSet.")
  }
  
  if (length(seqs) == 0L) {
    stop("Reference FASTA contains no sequences.")
  }
  
  # Remove Ns and gaps while preserving headers
  headers <- names(seqs)
  
  seqs <- Biostrings::DNAStringSet(
    gsub("[Nn.-]", "", as.character(seqs))
  )
  names(seqs) <- headers
  
  # Identify malformed records
  malformed <- (
    is.na(headers) |
      headers == "" |
      !grepl(";", headers, fixed = TRUE) |
      grepl("(^;|;;|;$)", headers) |
      grepl("(^|;)NA(;|$)", headers, ignore.case = TRUE) |
      Biostrings::width(seqs) == 0L |
      !grepl("[ACGT]", as.character(seqs), ignore.case = TRUE)
  )
  
  if (any(malformed)) {
    bad_headers <- headers[malformed]
    bad_headers[is.na(bad_headers) | bad_headers == ""] <- "<missing header>"
    
    report <- paste0(sum(malformed)," malformed FASTA record(s) detected:\n  ",
                     paste(head(bad_headers, 20L), collapse = "\n  "),
                     if (sum(malformed) > 20L) {
                       paste0(
                         "\n  ... and ",
                         sum(malformed) - 20L,
                         " additional record(s)"
                       )
                     } else {
                       ""
                     }
    )
    if (!drop_malformed) {
      stop(report, "\nSet 'drop_malformed = TRUE' to remove them.", call. = FALSE)
    }
    warning(report,"\nThese records will be removed.",call. = FALSE,immediate. = TRUE)
    seqs <- seqs[!malformed]
  }
  
  if (length(seqs) < 2L) {
    stop("Fewer than two usable reference sequences remain.")
  }
  
  if (orient) {
    if (!quiet) message("Orienting sequences")
    seqs <- DECIPHER::OrientNucleotides(seqs)
  }
  
  # Extract taxonomy and add Root rank
  taxonomy <- sub("^[^;]+;", "Root;", names(seqs))
  
  if (!quiet) {
    message(length(unique(taxonomy)), " unique terminal taxa")
  }
  
  # Prune over-represented terminal taxa
  group_counts <- table(taxonomy)
  remove <- logical(length(seqs))
  
  for (group in names(group_counts[group_counts > max_group_size])) {
    index <- which(taxonomy == group)
    remove[sample(index, length(index) - max_group_size)] <- TRUE
  }
  
  if (!quiet) {
    message(
      sum(remove),
      " sequences pruned from over-represented groups"
    )
  }
  
  train <- seqs[!remove]
  taxonomy <- taxonomy[!remove]
  groups <- taxonomy
  original_groups <- unique(groups)
  
  # Iteratively train and remove problem sequences
  remove <- logical(length(train))
  previous_problem_sequences <- integer()
  problem_sequences <- integer()
  
  for (i in seq_len(max_iterations)) {
    if (!quiet) message("Training iteration: ", i)
    
    training_set <- DECIPHER::LearnTaxa(
      train = train[!remove],
      taxonomy = taxonomy[!remove],
      rank = NULL
    )
    
    problem_sequences <- training_set$problemSequences$Index
    
    if (
      length(problem_sequences) == 0L ||
      identical(problem_sequences, previous_problem_sequences) ||
      i == max_iterations
    ) {
      break
    }
    
    previous_problem_sequences <- problem_sequences
    index <- which(!remove)[problem_sequences]
    remove[index] <- TRUE
    
    if (!allow_group_removal) {
      missing_groups <- setdiff(
        original_groups,
        groups[!remove]
      )
      
      restore <- index[groups[index] %in% missing_groups]
      remove[restore] <- FALSE
    }
    
    if (sum(!remove) < 2L) {
      stop("Fewer than two sequences remain after removing problematic training sequences.")
    }
  }
  
  if (!quiet) {
    message(sum(remove), " problem sequences removed")
    message(length(problem_sequences), " problem sequences remaining")
    message("Trained IDTAXA on ",sum(!remove)," sequences in ",format(Sys.time() - time, digits = 2))
  }
  
  return(training_set)
}


# Convert DNABin to DNAStringSet
DNAbin2DNAstringset <- function (x, remove_gaps = FALSE) {
  if(!methods::is(x, "DNAbin")){
    stop("Input must be a DNAbin")
  }
  x %>%
    as.list() %>%
    as.character %>%
    purrr::map(function(y){
      y[!y %in% tolower(Biostrings::DNA_ALPHABET)] <- "N"
      if(isTRUE(remove_gaps)){
        y[y=="-"] <- ""
      }
      paste0(y, collapse="")
    })%>%
    unlist %>%
    Biostrings::DNAStringSet()
}

# Convert DNA character string to binary DNAbin
char2DNAbin <- function (x) {
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80,
                     160, 112, 224, 176, 208, 240, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66,
               86, 72, 68, 78, 73, 45, 63)
  vec <- raw(89)
  vec[indices] <- dbytes
  s2d1 <- function(s) vec[as.integer(charToRaw(s))]
  out <- lapply(x, s2d1)
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  class(out) <- "DNAbin"
  return(out)
}

# Convert binary DNAbin to character
DNAbin2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77,
                     66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80, 160,
               112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if (is.list(x)) {
    out <- sapply(x, function(s) rawToChar(vec[as.integer(s)]))
  } else {
    out <- rawToChar(vec[as.integer(x)])
  }
  attr(out, "rerep.names") <- attr(x, "rerep.names")
  attr(out, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(out)
}

