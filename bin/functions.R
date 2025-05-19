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
    data.table::melt(id.vars = c("seq_name"), variable.name = "sample_id", 
                     value.name = "Abundance", variable.factor = FALSE)
  
  # Remove observations with no abundance
  df <- df[Abundance > 0]
  
  #sample data
  if(!is.null(phyloseq::sample_variables(ps))) {
    samdf <- phyloseq::sample_data(ps) %>%
      as("data.frame") %>% 
      data.table::as.data.table(keep.rownames = FALSE)
    df <- df[samdf, on = .(sample_id = sample_id)]
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
    step <- round(max(sample_sums(ps)) / 100)
  } else if (is.integer(step)){
    step <- step
  } else {
    stop("Step must be an integer or 'auto' ")
  }
  ps <- ps %>%
    phyloseq::prune_samples(sample_sums(.)>0, .) %>% 
    phyloseq::filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
  rare <- otu_table(ps) %>%
    as("matrix") %>%
    vegan::rarecurve(step=step) %>% 
    purrr::set_names(sample_names(ps)) %>%
    purrr::map_dfr(., function(x){
      b <- as.data.frame(x)
      b <- data.frame(OTU = b[,1], count = rownames(b))
      b$count <- as.numeric(gsub("N", "",  b$count))
      return(b)
    },.id="sample_id") %>%
    left_join(phyloseq::sample_data(ps)%>%
                as("matrix") %>%
                tibble::as_tibble() %>%
                dplyr::select(sample_id, fcid) %>%
                dplyr::distinct())
  
  gg.rare <- rare %>%
    ggplot2::ggplot() +
    geom_line(aes(x = count, y = OTU, group=sample_id), alpha=0.3)+
    geom_point(data = rare %>% 
                 group_by(sample_id) %>% 
                 top_n(1, count),
               aes(x = count, y = OTU, colour=count > threshold)) +
    scale_x_continuous(label = scales::label_number(scale_cut = append(scales::cut_short_scale(), 1, 1)))+
    scale_colour_manual(values=c("FALSE" = "#F8766D", "TRUE"="#619CFF"))+
    facet_wrap(fcid~., scales="free", ncol=1)+
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