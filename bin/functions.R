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
merge_seqtabs <- function(x) {
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
      stop("Failed to match all taxa from count matrix ", i, " to the merged taxon set.")
    }

    column_end <- column_start + ncol(current) - 1L
    merged[row_index, seq.int(column_start, column_end)] <- current
    column_start <- column_end + 1L
  }

  merged
}

merge_taxtabs <- function(x) {

    merged <- data.table::rbindlist(
        lapply(
            x,
            function(current) {
                data.table::as.data.table(
                    current,
                    keep.rownames = "seq_name"
                )
            }
        ),
        use.names = TRUE,
        fill = TRUE
    )

    merged <- merged[!duplicated(seq_name)]

    seq_names <- merged$seq_name
    merged[, seq_name := NULL]

    merged <- as.matrix(merged)
    rownames(merged) <- seq_names

    merged
}

merge_samdfs <- function(x) {

    merged <- data.table::rbindlist(
        lapply(
            x,
            function(current) {
                data.table::as.data.table(
                    current,
                    keep.rownames = "sample_primers"
                )
            }
        ),
        use.names = TRUE,
        fill = TRUE
    )

    sample_ids <- merged$sample_primers
    merged[, sample_primers := NULL]

    merged <- as.data.frame(merged)
    rownames(merged) <- sample_ids

    merged
}

# Merge phyloseq objects containing disjoint samples.
#
# All OTU tables are standardized to taxa as rows before merging.
# Taxonomy and reference sequences for duplicate ASVs retain the first
# occurrence, subject to the behaviour of merge_taxtabs() and merge_refseq().
merge_phyloseq_disjoint_samples <- function(arguments) {

    if (!length(arguments)) {
        stop("arguments must contain at least one phyloseq object.")
    }

    arguments <- unname(arguments)

    # Validate input types.
    valid_inputs <- vapply(
        arguments,
        methods::is,
        logical(1),
        class2 = "phyloseq"
    )

    if (!all(valid_inputs)) {
        stop(
            "All inputs must be phyloseq objects. ",
            "Invalid input indices: ",
            paste(which(!valid_inputs), collapse = ", ")
        )
    }

    # Require sample names to be unique across all input objects.
    samples <- unlist(
        lapply(
            arguments,
            phyloseq::sample_names
        ),
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

    # Extract OTU tables and standardize to:
    # rows = ASVs, columns = samples.
    seqtabs <- lapply(
        arguments,
        function(ps) {
            seqtab <- as(phyloseq::otu_table(ps),"matrix")
            if (!phyloseq::taxa_are_rows(ps)) {seqtab <- t(seqtab)}
            storage.mode(seqtab) <- "integer"
            seqtab
        }
    )

    # Extract taxonomy tables.
    taxtabs <- lapply(
        arguments,
        function(ps) {
            as(phyloseq::tax_table(ps),"matrix")
        }
    )

    # Extract sample metadata.
    samdfs <- lapply(
        arguments,
        function(ps) {
            as(phyloseq::sample_data(ps),"data.frame")
        }
    )

    # Extract reference sequences.
    refseqs <- lapply(
        arguments,
        function(ps){
          phyloseq::refseq(ps)
        }
    )

    # Merge the individual components.
    seqtab_merged <- merge_seqtabs(seqtabs)
    taxtab_merged <- merge_taxtabs(taxtabs)
    samdf_merged <- merge_samdf(samdfs)
    refseq_merged <- merge_refseq(refseqs)

    # Synchronize all taxa-bearing components with the merged abundance table.
    taxa <- rownames(seqtab_merged)

    missing_taxa <- list(
        taxtab = setdiff(taxa, rownames(taxtab_merged)),
        refseq = setdiff(taxa, names(refseq_merged))
    )

    missing_taxa <- missing_taxa[
        lengths(missing_taxa) > 0L
    ]

    if (length(missing_taxa)) {
        details <- vapply(
            names(missing_taxa),
            function(component) {
                paste0(
                    component,
                    ": ",
                    paste(missing_taxa[[component]], collapse = ", ")
                )
            },
            character(1)
        )

        stop(
            "Merged components are missing ASVs present in the sequence table. ",
            paste(details, collapse = "; ")
        )
    }

    taxtab_merged <- taxtab_merged[taxa, , drop = FALSE]
    refseq_merged <- refseq_merged[taxa]

    # Synchronize sample metadata with the merged abundance table.
    samples <- colnames(seqtab_merged)

    missing_samples <- setdiff(
        samples,
        rownames(samdf_merged)
    )

    if (length(missing_samples)) {
        stop(
            "Merged sample metadata is missing samples present in the ",
            "sequence table: ",
            paste(missing_samples, collapse = ", ")
        )
    }

    samdf_merged <- samdf_merged[samples, , drop = FALSE]

    # Construct the merged phyloseq object once.
    phyloseq::phyloseq(
        phyloseq::otu_table(seqtab_merged, taxa_are_rows = TRUE),
        phyloseq::tax_table(taxtab_merged),
        phyloseq::sample_data(samdf_merged),
        phyloseq::refseq(refseq_merged)
    )
}

summarise_phyloseq <- function(
    ps = NULL,
    otutab = NULL,
    taxtab = NULL,
    seqs = NULL
) {

    # Extract components from phyloseq if supplied.
    if (!is.null(ps)) {
        otutab <- as(phyloseq::otu_table(ps),"matrix")
        storage.mode(otutab) <- "integer"

        if (!phyloseq::taxa_are_rows(ps)) {
            otutab <- t(otutab)
        }

        taxtab <- data.table::as.data.table(
            as(phyloseq::tax_table(ps),"matrix"),
            keep.rownames = "seq_name"
        )

        seqs <- phyloseq::refseq(ps)

    } else {

        if (is.null(otutab)) {
            stop("Either ps or otutab must be supplied.")
        }

        if (!is.matrix(otutab)) {
            stop("otutab must be a matrix.")
        }
    }

    # Convert abundance matrix to wide output table.
    sample_cols <- colnames(otutab)

    df <- data.table::as.data.table(
        otutab,
        keep.rownames = "seq_name"
    )

    # Add sequences if available.
    if (!is.null(seqs)) {
        df[,sequence := as.character(seqs[seq_name])]
    }

    # Add taxonomy if available.
    if (!is.null(taxtab)) {
        taxtab <- data.table::as.data.table(
            taxtab,
            keep.rownames = "seq_name"
        )

        tax_cols <- setdiff(
            names(taxtab),
            "seq_name"
        )

        df <- taxtab[
            df,
            on = "seq_name"
        ]

        data.table::setcolorder(
            df,
            c(
                "seq_name",
                if (!is.null(seqs)) "sequence",
                tax_cols,
                sample_cols
            )
        )
    }

    df
}

# Convert a phyloseq object into a long data.table containing one row per non-zero sample–taxon observation.
melt_phyloseq <- function(
    ps = NULL,
    otutab = NULL,
    taxtab = NULL,
    samdf = NULL,
    seqs = NULL
) {

    # Extract components from phyloseq if supplied.
    if (!is.null(ps)) {
        otutab <- as(phyloseq::otu_table(ps),"matrix")
        storage.mode(otutab) <- "integer"

        if (!phyloseq::taxa_are_rows(ps)) {
            otutab <- t(otutab)
        }

        if (length(phyloseq::sample_variables(ps)) > 0L) {
            samdf <- data.table::as.data.table(
                as(phyloseq::sample_data(ps),"data.frame"),
                keep.rownames = "sample_primers"
            )
        }

        if (length(phyloseq::rank_names(ps)) > 0L) {
            taxtab <- data.table::as.data.table(
                as(phyloseq::tax_table(ps),"matrix"),
                keep.rownames = "seq_name"
            )
        }

        seqs <- phyloseq::refseq(ps)

    } else {
        # Individual objects input
        if (is.null(otutab)) {
            stop("Either ps or otutab must be supplied.")
        }

        if (!is.matrix(otutab)) {
          stop("otutab must be a matrix.")
        }
    }

    # Identify non-zero observations.
    index <- which(otutab > 0, arr.ind = TRUE )

    # Exit early if no observations exist.
    if (!nrow(index)) {
        return(data.table::data.table())
    }

    # Convert to one row per sample-taxon observation.
    df <- data.table::data.table(
        seq_name = rownames(otutab)[index[, "row"]],
        sample_primers = colnames(otutab)[index[, "col"]],
        Abundance = otutab[index]
    )

    # Add sample metadata.
    if (!is.null(samdf)) {
        samdf <- data.table::as.data.table(
            samdf,
            keep.rownames = "sample_primers"
        )

        sample_columns <- setdiff(
            names(samdf),
            "sample_primers"
        )

        if (length(sample_columns)) {
            sample_index <- match(
                df$sample_primers,
                samdf$sample_primers
            )

            df <- cbind(
                df,
                samdf[
                    sample_index,
                    sample_columns,
                    with = FALSE
                ]
            )
        }
    }

    # Add taxonomy.
    if (!is.null(taxtab)) {
        taxtab <- data.table::as.data.table(
            taxtab,
            keep.rownames = "seq_name"
        )

        taxonomy_columns <- setdiff(
            names(taxtab),
            "seq_name"
        )

        if (length(taxonomy_columns)) {
            taxon_index <- match(
                df$seq_name,
                taxtab$seq_name
            )

            df <- cbind(
                df,
                taxtab[
                    taxon_index,
                    taxonomy_columns,
                    with = FALSE
                ]
            )
        }
    }

    # Add sequences.
    if (!is.null(seqs)) {
        sequence_index <- match(df$seq_name,names(seqs))

        df[
            ,
            sequence := as.character(
                seqs[sequence_index]
            )
        ]
    }

    data.table::setorder(
        df,
        -Abundance,
        seq_name
    )

    df
}

# Rarefaction / species accumulation plot
rareplot <- function(
      ps = NULL,
      otutab = NULL,
      samdf = NULL,
      step = "auto",
      threshold = 0
  ) {

    # Extract components from phyloseq if supplied.
    if (!is.null(ps)) {
        otutab <- as(phyloseq::otu_table(ps), "matrix")
        storage.mode(otutab) <- "integer"

        if (!phyloseq::taxa_are_rows(ps)) {
            otutab <- t(otutab)
        }

        samdf <- as(phyloseq::sample_data(ps),"data.frame")

    } else {
        if (is.null(otutab)) {
            stop("Either ps or otutab must be supplied.")
        }
        if (!is.matrix(otutab)) {
            stop("otutab must be a matrix.")
        }
        if (is.null(samdf)) {
            stop("samdf must be supplied when ps is not provided.")
        }
    }

    # Remove empty taxa and samples.
    otutab <- otutab[
        rowSums(otutab) > 0,
        colSums(otutab) > 0,
        drop = FALSE
    ]

    # Determine rarefaction step size.
    if (identical(step, "auto")) {
        step <- max(1L, ceiling(median(colSums(otutab)) / 100))
    } else if (is.numeric(step)) {
        step <- as.numeric(ceiling(step))
    } else {
        stop("Step must be numeric or 'auto'.")
    }

    # vegan expects samples as rows and taxa as columns.
    rare <- vegan::rarecurve(
        t(otutab),
        step = step,
        tidy = TRUE
    ) %>%
        dplyr::rename(
            sample_primers = Site,
            count = Sample,
            ASV = Species
        ) %>%
        dplyr::left_join(
            samdf %>%
                dplyr::select(
                    sample_primers,
                    read_group
                ),
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


# Convert missing lower-rank classifications to labels based on the
# lowest successfully classified rank.
na_to_unclassified <- function(
    x,
    sep = "__",
    rownames = TRUE,
    quiet = FALSE
) {

  # Extract the taxonomy table and record the input type.
  if (inherits(x, "phyloseq")) {
    tax <- as(
      phyloseq::tax_table(x),
      "matrix"
    )

    input_type <- "phyloseq"
    id_col <- NULL

    if (!quiet) {
      message("Input is a phyloseq object")
    }

  } else if (is.matrix(x)) {
    tax <- x
    input_type <- "matrix"
    id_col <- NULL

  } else if (is.data.frame(x)) {
    input_type <- "data.frame"

    # Detect a supported identifier column.
    id_col <- intersect(
      c("OTU", "seq_name"),
      names(x)
    )[1L]

    if (!is.na(id_col)) {
      tax_names <- as.character(x[[id_col]])

      tax <- as.matrix(
        x[
          ,
          setdiff(names(x), id_col),
          drop = FALSE
        ]
      )

      rownames(tax) <- tax_names

    } else {
      id_col <- NULL
      tax <- as.matrix(x)
    }

  } else {
    stop(
      "x must be a matrix, data frame, or phyloseq object."
    )
  }

  # Ensure taxonomy values are stored as character strings.
  storage.mode(tax) <- "character"

  if (is.null(rownames(tax))) {
    stop(
      "The taxonomy table must have row names, ",
      "an OTU column, or a seq_name column."
    )
  }

  if (is.null(colnames(tax))) {
    stop(
      "The taxonomy table must have named taxonomy columns."
    )
  }

  # Treat empty strings as missing classifications.
  tax[tax == ""] <- NA_character_

  # Identify which taxonomy values are classified.
  classified <- !is.na(tax)

  # Identify the lowest classified rank in each row.
  lowest_index <- max.col(
    classified *
      rep(
        seq_len(ncol(tax)),
        each = nrow(tax)
      ),
    ties.method = "last"
  )

  # Identify rows containing at least one classification.
  has_classification <- rowSums(classified) > 0L
  row_index <- seq_len(nrow(tax))

  # Create one propagated label per row from its lowest classified rank.
  replacement <- rep(
    NA_character_,
    nrow(tax)
  )

  replacement[has_classification] <- paste0(
    substr(
      colnames(tax)[lowest_index[has_classification]],
      1L,
      1L
    ),
    sep,
    tax[
      cbind(
        row_index[has_classification],
        lowest_index[has_classification]
      )
    ]
  )

  # Locate missing taxonomy values in classified rows.
  missing <- which(
    is.na(tax) &
      has_classification[row(tax)],
    arr.ind = TRUE
  )

  # Fill missing ranks using the row-specific propagated label.
  if (nrow(missing)) {
    tax[missing] <- replacement[
      missing[, "row"]
    ]
  }

  # Return a modified phyloseq object.
  if (input_type == "phyloseq") {
    phyloseq::tax_table(x) <- phyloseq::tax_table(tax)
    return(x)
  }

  # Return a matrix when the input was a matrix and row names are requested.
  if (input_type == "matrix" && rownames) {
    return(tax)
  }

  # Convert the final taxonomy matrix to a data frame.
  tax_df <- as.data.frame(
    tax,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Return a data frame with taxonomy identifiers in the row names.
  if (rownames) {
    return(tax_df)
  }

  # Restore the detected identifier column, defaulting to OTU.
  output_id_col <- if (is.null(id_col)) {
    "OTU"
  } else {
    id_col
  }

  output <- data.frame(
    identifier = rownames(tax),
    tax_df,
    row.names = NULL,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  names(output)[1L] <- output_id_col

  output
}


# Convert labels generated by na_to_unclassified() back to NA.
unclassified_to_na <- function(
    x,
    sep = "__",
    rownames = TRUE,
    quiet = FALSE
) {

  # Extract the taxonomy table and record the input type.
  if (inherits(x, "phyloseq")) {
    tax <- as(phyloseq::tax_table(x), "matrix")
    input_type <- "phyloseq"

    if (!quiet) {
      message("Input is a phyloseq object")
    }

  } else if (is.matrix(x)) {
    tax <- x
    input_type <- "matrix"

  } else if (is.data.frame(x)) {
    input_type <- "data.frame"

    if ("OTU" %in% names(x)) {
      otu_names <- as.character(x[["OTU"]])
      tax <- as.matrix(x[, setdiff(names(x), "OTU"), drop = FALSE])
      rownames(tax) <- otu_names
    } else {
      tax <- as.matrix(x)
    }

  } else {
    stop("x must be a matrix, data frame, or phyloseq object.")
  }

  # Ensure taxonomy values are character rather than factors.
  storage.mode(tax) <- "character"

  if (is.null(rownames(tax))) {
    stop("The taxonomy table must have row names or an OTU column.")
  }

  # Match labels such as G__Drosophila or F__Drosophilidae.
  separator_pattern <- paste0(
    "^[[:alnum:]]",
    escape_regex(sep)
  )

  generated <- grepl(
    separator_pattern,
    tax
  )

  generated[is.na(generated)] <- FALSE
  tax[generated] <- NA_character_

  # Return the same broad object type as the input.
  if (input_type == "phyloseq") {
    phyloseq::tax_table(x) <- phyloseq::tax_table(tax)
    return(x)
  }

  if (input_type == "matrix") {
    if (rownames) {
      return(tax)
    }

    return(data.frame(
      OTU = rownames(tax),
      tax,
      row.names = NULL,
      check.names = FALSE
    ))
  }

  tax <- as.data.frame(
    tax,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (rownames) {
    return(tax)
  }

  data.frame(
    OTU = rownames(tax),
    tax,
    row.names = NULL,
    check.names = FALSE
  )
}


# Escape a literal separator before using it in a regular expression.
escape_regex <- function(x) {
  gsub(
    "([][{}()+*^$|\\\\?.])",
    "\\\\\\1",
    x
  )
}

# Validate inputs are suitable for creating a phyloseq object.
validate_phyloseq_inputs <- function(
    otutab,
    taxtab,
    samdf,
    seqs = NULL,
    tax_cols = NULL
) {

   # Report identifiers that differ between two objects.
    check_names <- function(x, y, x_name, y_name) {
        only_x <- setdiff(x, y)
        only_y <- setdiff(y, x)

        if (length(only_x) || length(only_y)) {
            stop(x_name, " and ", y_name," contain different identifiers.",
                if (length(only_x)) {
                    paste0("\nPresent only in ", x_name,": ",paste(only_x, collapse = ", "))
                },
                if (length(only_y)) {
                    paste0("\nPresent only in ", y_name, ": ",paste(only_y, collapse = ", "))
                }
            )
        }
    }

    # Basic object checks
    if (!is.matrix(otutab)) {stop("otutab must be a matrix.")}
    if (!is.matrix(taxtab)) {stop("taxtab must be a matrix.")}
    if (!is.data.frame(samdf)) {stop("samdf must be a data.frame.")}
    if (is.null(rownames(otutab))) {stop("otutab must have taxa as rownames.")}
    if (is.null(colnames(otutab))) {stop("otutab must have sample names as colnames.")}
    if (is.null(rownames(taxtab))) {stop("taxtab must have taxa as rownames.")}
    if (is.null(rownames(samdf))) {stop("samdf must have sample names as rownames.")}
    if (anyDuplicated(rownames(samdf))) {stop("samdf contains duplicated sample names.")}

    # Check taxonomy columns.

    if (!is.null(tax_cols) && !identical(colnames(taxtab), tax_cols)) {
        stop( "Unexpected taxonomy columns.\n",
            "Expected: ",paste(tax_cols, collapse = ", "),
            "\nObserved: ",paste(colnames(taxtab), collapse = ", ")
        )
    }

    # Check sample names between the OTU table and sample metadata.
    check_names(
        colnames(otutab),
        rownames(samdf),
        "otutab samples",
        "samdf samples"
    )

    # Check taxon names between the OTU and taxonomy tables.
    check_names(
        rownames(otutab),
        rownames(taxtab),
        "otutab taxa",
        "taxtab taxa"
    )

    # Check taxon names between the OTU table and reference sequences.
    if (!is.null(seqs)) {
        if (is.null(names(seqs))) {
            stop("seqs must have taxon names.")
        }

        check_names(
            rownames(otutab),
            names(seqs),
            "otutab taxa",
            "seqs"
        )
    }

    invisible(TRUE)
}