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

# Sample validation -------------------------------------------------------


#Update the sample sheet and logging sheet to deal with any newly demultiplexed files
step_demux_samdf <- function(samdf){
    out <- samdf %>%
      dplyr::group_by(sample_id) %>%
      dyplr::group_split() %>%
      purrr::map(function(x){
        if(any(stringr::str_detect(x$pcr_primers, ";"), na.rm = TRUE)){
          primer_names <- unlist(stringr::str_split(unique(x$pcr_primers), ";")) 
          x <- x %>% 
            dplyr::mutate(count = length(primer_names)) %>% #Replicate the samples
            tidyr::uncount(count) %>%
            dplyr::mutate(pcr_primers = unlist(stringr::str_split(unique(x$pcr_primers), ";")),
                   for_primer_seq = unlist(stringr::str_split(unique(x$for_primer_seq), ";")),
                   rev_primer_seq = unlist(stringr::str_split(unique(x$rev_primer_seq), ";")),
                   sample_id = paste0(sample_id, "_",pcr_primers)
            ) 
        }
        if(!all(stringr::str_detect(x$sample_id ,paste0(x$pcr_primers, "$")))){
          x <- x %>%
            dplyr::mutate(sample_id = paste0(sample_id, "_",pcr_primers))
        }
        return(x)
      }) %>%
      dplyr::bind_rows()
    
    # Check if files exist
    data_folders <- paste0(list.dirs("data", recursive=FALSE), "/01_trimmed")
    fastqFs <- purrr::map(data_folders,list.files, pattern="_R1_", full.names = TRUE) %>%
      unlist() %>%
      stringr::str_remove(pattern = "^(.*)\\/") %>%
      stringr::str_remove(pattern = "(?:.(?!_S))+$")
    fastqFs <- fastqFs[!stringr::str_detect(fastqFs, "Undetermined")]
    #Check missing fastqs
    if (length(setdiff(out$sample_id, fastqFs)) > 0) {
      warning(paste0("The fastq file: ",
                     setdiff(out$sample_id, fastqFs),
                     " is missing, dropping from samplesheet \n")) 
      out <- out %>%
        filter(!sample_id %in% setdiff(out$sample_id, fastqFs))
    }
    return(out)
}


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


step_phyloseq <- function(seqtab, taxtab, samdf, seqs=NULL, phylo=NULL, name_variants=FALSE){

  # Check if seqtab is a path
  if(is(seqtab, "character")){
    if(file.exists(seqtab)){
      seqtab <- readRDS(normalizePath(seqtab))
    } else {
      stop("seqtab does not exist")
    }
  }
  
  # Check if taxtab is a path
  if(is(taxtab, "character") ){
    if(file.exists(taxtab)){
      taxtab <- readRDS(normalizePath(taxtab))
    } else {
      stop("taxtab does not exist")
    }
  }
  
  # Check if samdf is a path - if so read in
  if(is(samdf, "character")){
    if (file.exists(samdf)){
      samdf <- readr::read_csv(normalizePath(samdf))
    } else {
      stop("samdf does not exist")
    }
  } 
  # Check if samdf is a path - if so read in
  if(is(seqs, "character")){
    if (file.exists(seqs)){
      seqs <- readRDS(normalizePath(seqs))
      names(seqs) <- seqs
    } else {
      stop("seqs path does not exist")
    }
  } else if (class(seqs) == "DNAStringSet"){
    seqs <- seqs
  } else {
    seqs <- Biostrings::DNAStringSet(colnames(seqtab))
    names(seqs) <- seqs
  }
  
  #Check if phy is a path - if so read in
  if(is(phylo, "character")){
    if (file.exists(phylo)){
      phy <- ape::read.tree(normalizePath(phylo))
    } else {
      stop("phy path does not exist")
    }
  } else if (class(phylo) == "phylo"){
    phy <- phylo
  } else {
    phy <- NULL
  }
  
  #Extract start of sequence names
  rownames(seqtab) <- stringr::str_remove(rownames(seqtab), pattern="_S[0-9]+_R[1-2]_.*$")
  
  #Load sample information
  samdf <- samdf %>%
    dplyr::filter(!duplicated(sample_id)) %>%
    as.data.frame()%>%
    magrittr::set_rownames(.$sample_id)
  
  missing_seqtab_asvs <- length(colnames(seqtab)[!colnames(seqtab) %in% rownames(taxtab)])
  missing_taxtab_asvs <-length(rownames(taxtab)[!rownames(taxtab) %in% colnames(seqtab)])
  
  if(missing_seqtab_asvs > 0 | missing_taxtab_asvs){
    stop(paste0(missing_seqtab_asvs, " ASVs are in seqtab and not taxtab, and ",
                missing_taxtab_asvs, " ASVs are in taxtab but not seqtab"))
  }
  
  if(is.null(phy)){
    ps <- phyloseq(phyloseq::tax_table(taxtab),
                   phyloseq::sample_data(samdf),
                   phyloseq::otu_table(seqtab, taxa_are_rows = FALSE),
                   phyloseq::refseq(seqs))
  } else {
    ps <- phyloseq(phyloseq::tax_table(taxtab),
                   phyloseq::sample_data(samdf),
                   phyloseq::otu_table(seqtab, taxa_are_rows = FALSE),
                   phy_tree(phy),
                   phyloseq::refseq(seqs))
  }
  
  if(name_variants){
    phyloseq::taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)),"-",tax_table(ps)[,8])
  }
  
  if(nrow(seqtab) > nrow(phyloseq::sample_data(ps))){
    message("Warning: the following samples were not included in phyloseq object, check sample names match the sample metadata")
    message(rownames(seqtab)[!rownames(seqtab) %in% phyloseq::sample_names(ps)])
  }
  return(ps)
}

new_bind <- function(a, b) {
  common_cols <- intersect(names(a), names(b))
  b[common_cols] <- map2_df(b[common_cols], 
                            map(a[common_cols], class), ~{class(.x) <- .y;.x})
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
merge_phyloseq_new <- function (arguments){
  comp.list <- list()
  for (i in 1:length(arguments)) {
    comp.list <- c(comp.list, phyloseq:::splat.phyloseq.objects(arguments[[i]]))
  }
  merged.list <- list()
  for (i in unique(names(comp.list))) {
    i.list <- comp.list[names(comp.list) == i]
    if (length(i.list) == 1) {
      merged.list <- c(merged.list, i.list)
    } else {
      x1 <- i.list[[1]]
      for (j in 2:length(i.list)) {
        x1 <- phyloseq::merge_phyloseq_pair(x1, i.list[[j]])
      }
      x1 <- list(x1)
      names(x1) <- i
      merged.list <- c(merged.list, x1)
    }
    gc()
  }
  names(merged.list) <- NULL
  return(do.call(phyloseq, merged.list))
}

summarise_phyloseq <- function(ps){
  phyloseq::otu_table(ps) %>%
    # t() %>%     # removed because rows are sequence names not samples
    as.data.frame() %>%
    tibble::rownames_to_column("seq_name") %>%
    dplyr::left_join(
      phyloseq::tax_table(ps) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("seq_name"),
      by = "seq_name"
    ) %>%  
    dplyr::left_join(
      phyloseq::refseq(ps) %>% 
        as.character() %>% 
        tibble::enframe(name="seq_name", value="sequence"),
      by = "seq_name"
    )  %>%
    dplyr::select(seq_name, sequence, phyloseq::rank_names(ps), phyloseq::sample_names(ps))
}

melt_phyloseq <- function(ps) {
  
  # Enforce OTU table orientation. Redundant-looking step
  seqtab <- phyloseq::otu_table(ps)
  if(!taxa_are_rows(seqtab)){seqtab <- t(seqtab)}
  
  # Convert otu table to tall form (one sample-taxon per row)
  df <- seqtab %>% 
    as("matrix") %>%
    data.table::as.data.table(keep.rownames = "seq_name") %>%
    data.table::melt(id.vars = c("seq_name"), variable.name = "sample_primers", 
                     value.name = "Abundance", variable.factor = FALSE)
  
  # Remove observations with no abundance
  df <- df[Abundance > 0]
  
  #sample data
  if(!is.null(phyloseq::sample_variables(ps))) {
    samdf <- phyloseq::sample_data(ps) %>%
      as("data.frame") %>% 
      data.table::as.data.table(keep.rownames = FALSE)
    df <- df[samdf, on = .(sample_primers = sample_primers)]
  }
  
  # Add the tax table if it exists
  if(!is.null(phyloseq::rank_names(ps))) {
    taxtab <- phyloseq::tax_table(ps) %>%
      as("matrix") %>%
      data.table::as.data.table(keep.rownames = "seq_name")
    df <- df[taxtab, on = .(seq_name = seq_name)]
  }
  
  # Add the sequences if they exist
  if(!is.null(phyloseq::refseq(ps))) {
    seqs <- phyloseq::refseq(ps) %>%
      as("character") %>%
      data.table::as.data.table(keep.rownames = "seq_name")
    data.table::setnames(seqs, old = c("seq_name", "."), new = c("seq_name", "sequence"))
    df <- df[seqs, on = .(seq_name = seq_name)]
  }
  
  # Arrange by Abundance, then sequence names (to approx. phyloseq behavior)
  df <- df %>%
    data.table::setorder(-Abundance, seq_name) 
  return(df)
}

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
  } else {
    otutab <- as(phyloseq::otu_table(ps), "matrix")
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

#' write fasta
#'
#' @param x a list of sequences in DNAbin or AAbin format, or a vector of sequences as concatenated upper-case character strings.
#' @param file character string giving a valid file path to output the text to. If file = "" (default setting) the text file is written to the console.
#' @param compress logical indicating whether the output file should be gzipped.
#' @param quiet Whether progress should be printed to consoe
#'
#' @return
#' @export
#'
#' @examples
write_fasta <- function(x, file = "", compress = FALSE, quiet=FALSE) {
  if(stringr::str_detect(file, "\\.gz$") & !compress){
    compress <- TRUE
    if(!quiet) message(".gz detected in filename, compressing output file")
  }
  if (!is.null(dim(x))) {
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if (inherits(x, "DNAbin")) {
    tmp <- DNAbin2char(x)
  } else if (is.list(x)) {
    if (length(x[[1]] == 1)) {
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
  if(!quiet){message("Wrote ", length(tmp), " sequences to ", file)}
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

#' Train IDTAXA
#'
#' @param seqs A DNAbin
#' @param max_group_size The maximum size of any taxonomic group. This can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
#' @param max_iterations The number of iterations to conduct training for in order to identify  any training sequences whose assigned classifications completely disagree with their predicted classification.
#' @param allow_group_removal  Whether sequences that are the last remaining representatives of an entire group in the training data can be removed.
#' This can occur if the entire group appears to be misplaced in the taxonomic tree.
#' @param orient Training sequences must all be in the same orientation. Set this to TRUE to reorient the sequences if you are unsure.
#' @param quiet Whether progress should be printed to console.
#'
#' @return
#' @export
#' @import stringr
#' @importFrom DECIPHER RemoveGaps
#' @importFrom DECIPHER OrientNucleotides
#' @importFrom DECIPHER LearnTaxa
#'
train_idtaxa <- function(seqs, max_group_size=10, max_iterations = 3,  allow_group_removal = TRUE, orient=FALSE, quiet = FALSE) {
  time <- Sys.time() # get time

  #Remove NA's
  if(length(names(seqs)[stringr::str_detect(names(seqs), ";NA;")]) > 1){
    remove <- names(seqs)[stringr::str_detect(names(seqs), ";NA;")]
    subset <- seqs[!names(seqs) %in% remove]
    message(paste0(length(seqs) - length(subset)," Sequences with NA's in taxonomy removed"))
    seqs <- subset
  }

  # Convert to DNAstringset for DECIPHER
  if(methods::is(seqs, "DNAbin")){
    seqs <- DNAbin2DNAstringset(seqs)
  }

  # Remove gaps
  if(!quiet){ message("Removing gaps")}
  seqs <- DECIPHER::RemoveGaps(seqs)

  if(orient){
    if(!quiet){ message("Orienting sequences")}
    seqs <- DECIPHER::OrientNucleotides(seqs)
  }

  # As taxonomies are encoded in the sequence names rather than a separate file, use:
  taxid <- NULL

  #Add Root Rank
  names(seqs) <- names(seqs)  %>% stringr::str_replace(";[;]*", ";Root;")

  # obtain the taxonomic assignments
  groups <- names(seqs) # sequence names

  # assume the taxonomy begins with 'Root;'
  groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
  groupCounts <- table(groups)
  u_groups <- names(groupCounts) # unique groups

  if (!quiet) {message(paste0(length(u_groups), " unique species"))}

  # Pruning training set
  remove <- logical(length(seqs))
  for (i in which(groupCounts > max_group_size)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
                   max_group_size)
    remove[index[-keep]] <- TRUE
  }
  if (!quiet) {message(paste0(sum(remove), " sequences pruned from over-represented groups"))}

  # Training the classifier
  probSeqsPrev <- integer() # suspected problem sequences from prior iteration

  train <- seqs[!remove]
  taxonomy <- gsub("(.*)(Root;)", "\\2", names(train))
  remove <- logical(length(train))

  for (i in seq_len(max_iterations)) {
    if (!quiet) {message("Training iteration: ", i, "\n", sep="")}
    # train the classifier
    trainingSet <- DECIPHER::LearnTaxa(train[!remove], taxonomy[!remove], taxid)

    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
      if (!quiet) {message("No problem sequences remaining.\n")}
      break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
      if (!quiet) {message("Iterations converged.\n")}
      break
    }
    if (i==max_iterations)
      break
    probSeqsPrev <- probSeqs

    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allow_group_removal) {
      # replace any removed groups
      missing <- !(u_groups %in% groups[!remove])
      missing <- u_groups[missing]
      if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
      }
    }
  }
  if (!quiet) {message(paste0(sum(remove), " sequences removed"))}
  if (!quiet) {message(paste0(length(probSeqs), " problem sequences remaining"))}

  time <- Sys.time() - time
  if (!quiet) (message(paste0("trained IDTAXA on ", length(seqs), " sequences in ", format(time, digits = 2))))
  return(trainingSet)
}


#' Convert DNABin to DNAStringSet
#'
#' @param x a DNABin object
#' @param remove_gaps Whether gaps should be removed
#'
#' @return
#' @export
#' @import purrr
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods is
#'
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

#' DNAbin to char
#'
#' @param x A DNAbin object
#'
#' @return
#' @export
#'
#' @examples
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

# Multithread -------------------------------------------------------------
#' Setup multithreading
#'
#' @param multithread Number of cores
#' @param quiet Whether progress should be printed to console
#'
#' @return
#' @import future
#'
#' @examples
setup_multithread <- function(multithread, quiet=FALSE){
  ncores <- future::availableCores()
  if(isTRUE(multithread)){
    cores <- ncores-1
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multisession, workers=cores)
  } else if (is.numeric(multithread) & multithread > 1){
    cores <- multithread
    if(cores > ncores){
      cores <- ncores
      warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
    }
    if(!quiet){message("Multithreading with ", cores, " cores")}
    future::plan(future::multisession, workers=cores)
  } else if(isFALSE(multithread) | multithread==1){
    future::plan(future::sequential)
  } else (
    stop("Multithread must be a logical or numeric vector of the numbers of cores to use")
  )
}
