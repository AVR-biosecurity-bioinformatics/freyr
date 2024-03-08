#### NOTE:
#### The following variables are defined in the module script blocks:
#### projectDir : the project directory absolute path; location of main.nf and "working directory"
#### data_loc : the name of the data directory (usually "data" but could be "test_data")


# Sample validation -------------------------------------------------------


#Update the sample sheet and logging sheet to deal with any newly demultiplexed files
step_demux_samdf <- function(samdf){
    out <- samdf %>%
      dplyr::group_by(sample_id) %>%
      group_split() %>%
      purrr::map(function(x){
        if(any(stringr::str_detect(x$pcr_primers, ";"), na.rm = TRUE)){
          primer_names <- unlist(stringr::str_split(unique(x$pcr_primers), ";")) 
          x <- x %>% 
            dplyr::mutate(count = length(primer_names)) %>% #Replicate the samples
            uncount(count) %>%
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

step_add_params <- function(samdf, params){
  out <- samdf %>%
    seqateurs::coalesce_join(params, by="pcr_primers")
  return(out)
}

step_check_files <- function(samdf, files, col_name=NULL){
  # Get file names
  fastqFs <- files[stringr::str_detect(files, "_R1_")]
  fastqRs <- files[stringr::str_detect(files, "_R2_")]
  
  # Get file name to check
  namecheck <- basename(fastqFs) %>%
    stringr::str_remove(pattern = "^(.*)\\/") %>%
    stringr::str_remove(pattern = "(?:.(?!_S))+$")
  namecheck <- namecheck[!stringr::str_detect(namecheck, "Undetermined")]
  
  #Check missing in samplesheet
  if (length(setdiff(namecheck, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(namecheck, samdf$sample_id), " are not in the sample sheet") }
  
  #Check missing fastqs
  if (length(setdiff(samdf$sample_id, namecheck)) > 0) {
    warning(paste0("The fastq file: ",
                   setdiff(samdf$sample_id, namecheck),
                   " is missing, dropping from samplesheet \n")) 
    samdf <- samdf %>%
      filter(!sample_id %in% setdiff(samdf$sample_id, namecheck))
  }
  
  # Hash the file to make sure they havent changed
  fwd_col_name = paste0(col_name, "_fwd")
  rev_col_name = paste0(col_name, "_rev")
  out <- samdf %>%
    mutate(!!fwd_col_name := purrr::map_chr(sample_id,~{
      str_to_check <- basename(fastqFs) %>%
        stringr::str_remove(pattern = "^(.*)\\/") %>%
        stringr::str_remove(pattern = "(?:.(?!_S))+$")
      rlang::hash_file(fastqFs[str_to_check==.x])
    }))%>%
    mutate(!!rev_col_name := purrr::map_chr(sample_id,~{
      str_to_check <- basename(fastqRs) %>%
        stringr::str_remove(pattern = "^(.*)\\/") %>%
        stringr::str_remove(pattern = "(?:.(?!_S))+$")
      rlang::hash_file(fastqRs[str_to_check==.x])
    }))
  return(out)
}

step_validate_folders <- function(project_dir){
  # Required directories
  list("data",
       "reference",
       "output/logs",
       "output/results/unfiltered",
       "output/results/filtered",
       "output/rds",
       "output/temp",
       "sample_data") %>%
    purrr::map(normalizePath) %>%
    purrr::walk(function(x){
      if(!dir.exists(x)){
        dir.create(x, recursive=TRUE)
      }
    })
}

step_validate_samdf <- function(samdf, data_dir){
  # Check if sampleids contain fcid, if not; attatch
  samdf <- samdf %>%
    dplyr::mutate(sample_id = case_when(
      !stringr::str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
      TRUE ~ sample_id
    ))
  
  # Check if forward read files match samplesheet
  data_dir <- normalizePath(data_dir)
  fastqFs <- purrr::map(list.dirs(data_dir, recursive=FALSE),
                        list.files, pattern="_R1_", full.names = TRUE) %>%
    unlist() %>%
    stringr::str_remove(pattern = "^(.*)\\/") %>%
    stringr::str_remove(pattern = "(?:.(?!_S))+$")
  fastqFs <- fastqFs[!stringr::str_detect(fastqFs, "Undetermined")]
  #Check missing in samplesheet
  if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") }
  
  # Check if reverse files match samplesheet
  data_dir <- normalizePath(data_dir)
  fastqRs <- purrr::map(list.dirs(data_dir, recursive=FALSE),
                        list.files, pattern="_R2_", full.names = TRUE) %>%
    unlist() %>%
    stringr::str_remove(pattern = "^(.*)\\/") %>%
    stringr::str_remove(pattern = "(?:.(?!_S))+$")
  fastqRs <- fastqRs[!stringr::str_detect(fastqRs, "Undetermined")]
  
  #Check missing in samplesheet
  if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") }
  
  #Check missing fastqs
  if (length(setdiff(samdf$sample_id, fastqFs)) > 0) {
    warning(paste0("The fastq file: ",
                   setdiff(samdf$sample_id, fastqFs),
                   " is missing, dropping from samplesheet \n")) 
    samdf <- samdf %>%
      filter(!sample_id %in% setdiff(samdf$sample_id, fastqFs))
  }
  return(samdf)
}
  



# Quality control ---------------------------------------------------------

step_seq_qc <- function(fcid, quiet=FALSE, write_all=FALSE){
  
  seq_dir <- paste0(projectDir,"/",data_loc,"/",fcid) # note: uses variables defined in module script
  qc_dir <- paste0(projectDir,"/output/logs/",fcid)
  
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  if(!dir.exists(paste0(seq_dir,"/InterOp"))){
    warning("InterOp folder must be present to run quality checks")
    out <- tibble(fcid = fcid,
                  reads_pf = NA_integer_,
                  reads_total = NA_integer_)
    return(out)
  }
  if(!file.exists(paste0(seq_dir,"/RunInfo.xml"))){
    warning("RunInfo.xml must be present to run quality checks")
    out <- tibble(fcid = fcid,
                  reads_pf = NA_integer_,
                  reads_total = NA_integer_)
    return(out)
  }
  
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  ## Sequencing run quality check using savR
  fc <- savR::savR(seq_dir)
  
  # Ensure indices are present
  if(length(fc@parsedData)==0){
    warning(paste0("Flow cell metrics could not be parsed for", fcid, " skipping seq run qc"))
    out <- tibble(fcid = fcid,
                  reads_pf = NA_integer_,
                  reads_total = NA_integer_)
    return(out)
  }
  
  if(write_all){
    readr::write_csv(savR::correctedIntensities(fc), normalizePath(paste0(qc_dir, "/correctedIntensities.csv")))
    readr::write_csv(savR::errorMetrics(fc), normalizePath(paste0(qc_dir, "/errorMetrics.csv")))
    readr::write_csv(savR::extractionMetrics(fc), normalizePath(paste0(qc_dir, "/extractionMetrics.csv")))
    readr::write_csv(savR::qualityMetrics(fc), normalizePath(paste0(qc_dir, "/qualityMetrics.csv")))
    readr::write_csv(savR::tileMetrics(fc), normalizePath(paste0(qc_dir, "/tileMetrics.csv")))
  }
  
  gg.avg_intensity <- fc@parsedData[["savCorrectedIntensityFormat"]]@data %>%
    dplyr::group_by(tile, lane) %>%
    dplyr::summarise(Average_intensity = mean(avg_intensity), .groups="drop") %>% 
    dplyr::mutate(side = case_when(
      stringr::str_detect(tile, "^11") ~ "Top",
      stringr::str_detect(tile, "^21") ~ "Bottom"
    ))%>%
    ggplot2::ggplot(aes(x=lane, y=as.factor(tile), fill=Average_intensity)) +
    geom_tile() +
    facet_wrap(~side, scales="free") +
    scale_fill_viridis_c()
  
  pdf(file=paste0(fcid,"_flowcell_qc.pdf"), width = 11, height = 8 , paper="a4r")
    plot(gg.avg_intensity)
    savR::pfBoxplot(fc)
    for (lane in 1:fc@layout@lanecount) {
    savR::qualityHeatmap(fc, lane, 1:fc@directions)
    }
  try(dev.off(), silent=TRUE)

  if(!quiet){message("Flow cell quality metrics written to: ", qc_dir)}
  
  # Return the total number of reads and passing filter
  out <- fc@parsedData[["savTileFormat"]]@data %>%
    dplyr::filter(code %in% c(100,101)) %>%
    dplyr::mutate(code = case_when(
      code == 100 ~ "reads_total",
      code == 101 ~ "reads_pf"
    ),fcid = stringr::str_remove(fc@flowcell, "^.*-")) %>% 
    dplyr::group_by(fcid, code) %>%
    dplyr::summarise(reads = sum(value), .groups="drop") %>%
    tidyr::pivot_wider(names_from = code,
                       values_from = reads)
  return(out)
}

step_sample_qc <- function(sample_id, fcid, multithread=FALSE, quiet=FALSE){
  
  seq_dir <- paste0("data/", fcid, "/")
  qc_dir <- paste0("output/logs/", fcid,"/FASTQC" )
  
  print(qc_dir)
  print(seq_dir)
  
  fq <- paste0(seq_dir, sample_id, "*.fastq.gz")
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
  }
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  # Handle multithreading
  cores <- setup_multithread(multithread)
  
  # Define fastqc function
  fastqc <- function (fq, qc.dir = NULL, threads = 1, fastqc.path = "bin/FastQC") {
    if (is.null(qc.dir)) {
      qc.dir <- file.path(dirname(fq), "FASTQC")
    } 
    if(!dir.exists(qc.dir)) {dir.create(qc.dir, recursive = TRUE)}
    
    if (.Platform$OS.type == "unix") {
      fastqc.path <- file.path(fastqc.path, "fastqc")
      cmd <- paste0(fastqc.path, " ", fq, "  --threads ", 
                    threads, " --outdir ", qc.dir)
      result <- system(cmd)
    } else {
      #.threads <- paste0("--threads=", threads)
      #.qc.dir <- paste0("--outdir=", qc.dir)
      args <- paste(" -Xmx250m -cp", paste0(fastqc.path,";",fastqc.path,"/sam-1.103.jar;",fastqc.path,"/jbzip2-0.9.jar uk.ac.babraham.FastQC.FastQCApplication"),
                       fq, collapse=" ")
      result <- processx::run(command = "java", 
                              args = args, echo = TRUE, echo_cmd = TRUE, 
                              spinner = TRUE, windows_verbatim_args = TRUE, error_on_status = FALSE, 
                              cleanup_tree = TRUE)  
      outfiles <- c(fs::dir_ls(dirname(fq), glob = paste0("*",sample_id, "*.html")),
                    fs::dir_ls(dirname(fq), glob = paste0("*",sample_id, "*.zip")))
      fs::file_move(outfiles, qc.dir)
      }
    return(result)
  }
  
  # Run fastqc
  qc_out <- fastqc(fq = fq, qc.dir= qc_dir, fastqc.path = "bin/FastQC", threads=cores)

  # Check if errored
  if(stringr::str_detect(qc_out$stderr, "Error:")){
    stop(qc_out$stderr %>% str_remove("Error: "))
  }
}

step_multiqc <- function(fcid, quiet=FALSE){
  qc_dir <- paste0("output/logs/", fcid,"/FASTQC/" )
  # Check that required files exist
  if(!dir.exists(qc_dir)) {
    stop("qc_dir doesnt exist, check that the correct path was provided")
  }
  # Write out fastqc report
  if(!quiet){
    ngsReports::writeHtmlReport(qc_dir, overwrite = TRUE, gcType ="Genome",  quiet=quiet)
    message("Sample quality metrics written to: ", qc_dir)
  } else {
    suppressMessages(ngsReports::writeHtmlReport(qc_dir, overwrite = TRUE, gcType ="Genome",  quiet=quiet))
  }
  
}

step_switching_calc <- function(fcid, barcode_mismatch=1, quiet=FALSE){
  
  seq_dir <- paste0(projectDir,"/",data_loc,"/",fcid) # note: uses variables defined in module script
  qc_dir <- paste0(projectDir,"/output/logs/",fcid)
  
  # Check that required files exist
  if(!dir.exists(seq_dir)) {
    stop("input directory doesnt exist, check that the correct path was provided")
  }
  
  # Check if undetermined reads file exists
  if(!any(stringr::str_detect(list.files(seq_dir, pattern="_R1_", full.names = TRUE), "Undetermined"), na.rm = TRUE)){
    warning("Error, an Undetermined reads fastq must be present to calculate index switching")
    res <- tibble(expected = NA_integer_, observed=NA_integer_, switch_rate=NA_integer_)
    return(res)
  }
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  # Summarise indices
  indices <- sort(list.files(seq_dir, pattern="_R1_", full.names = TRUE)) 
  indices <- indices[stringr::str_detect(indices, ".fastq.gz$")] %>%
    purrr::set_names() %>%
    purrr::map(seqateurs::summarise_index) %>%
    dplyr::bind_rows(.id="Sample_Name")%>%
    dplyr::arrange(desc(Freq)) %>% 
    dplyr::mutate(Sample_Name = Sample_Name %>% 
                    stringr::str_remove(pattern = "^(.*)\\/") %>%
                    stringr::str_remove(pattern = "(?:.(?!_S))+$"))
  
  # Ensure indices are present
  if(all(is.na(indices$Freq))){
    warning(paste0("No index sequences present in fastq headers for run", fcid, " no switch rate calculated"))
    res <- tibble(expected = NA_integer_, observed=NA_integer_, switch_rate=NA_integer_)
    return(res)
  }
  
  combos <- indices %>% 
    dplyr::filter(!stringr::str_detect(Sample_Name, "Undetermined")) %>%
    dplyr::select(index, index2) %>%
    tidyr::expand(index, index2)
  
  #get unused combinations resulting from index switching
  switched <- combos %>%
    dplyr::left_join(indices, by=c("index", "index2")) %>%
    tidyr::drop_na()
  
  # Get a list of orignally applied indexes - Could get this from sample sheet instead
  applied_indices <- switched %>%
    dplyr::filter(!stringr::str_detect(Sample_Name, "Undetermined")) %>%
    dplyr::group_by(Sample_Name) %>%
    group_modify(~{
      .x %>%
        dplyr::top_n(n=1, Freq) %>%
        dplyr::slice(1) %>%  # Handle ties
        dplyr::mutate(Freq = sum(.x$Freq)) 
    })
  
  # Check if indices are combinatorial
  if(any(duplicated(applied_indices$index)) | any(duplicated(applied_indices$index2))){
    warning(paste0("Combinatorial indexes detected for", fcid, " no switch rate calculated"))
    res <- tibble(expected = NA_integer_, observed=NA_integer_, switch_rate=NA_integer_, contam_rate=NA_integer_)
    return(res)
  }
  
  # Get other undetermined reads which had completely unapplied indexes
  other_reads <- anti_join(indices,combos, by=c("index", "index2")) %>%
    dplyr::summarise(sum = sum(Freq, na.rm = TRUE)) %>%
    dplyr::pull(sum)
  
  #Summary of index switching rate
  if(any(stringr::str_detect(switched$Sample_Name, "Undetermined"), na.rm = TRUE)){
    res <- switched %>%
      dplyr::mutate(type = case_when(
        !stringr::str_detect(Sample_Name, "Undetermined") ~ "expected",
        stringr::str_detect(Sample_Name, "Undetermined") ~ "observed"
      )) %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(switch_rate = sum(Freq), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = type,
                         values_from = switch_rate) %>%
      dplyr::mutate(switch_rate =  observed / expected ) %>%
      dplyr::mutate(contam_rate =  switch_rate^2 )
  } else {
    res <- tibble(expected = NA_integer_, observed=NA_integer_, switch_rate=NA_integer_, contam_rate=NA_integer_)
    return(res)
  }
  
  if(!quiet){message("Index switching rate calculated as: ", res$switch_rate)}
  
  #Plot switching - handling barcode mismatch during demultiplexing 
 
  # Update indexes using their hamming distance to those originally appliex
  switch_plot_dat <- switched %>%
    dplyr::mutate(index = purrr::map(index, ~{
      index_list <- applied_indices$index
      index_dist <- stringdist::stringdist(.x,index_list)
      # Remove those above barcode_mismatch threshold
      index_list <- index_list[index_dist <= barcode_mismatch]
      index_dist <- index_dist[index_dist <= barcode_mismatch]
      return(index_list[which.min(index_dist)])
    })) %>%
    dplyr::mutate(index2 = purrr::map(index2, ~{
      index_list <- applied_indices$index2
      index_dist <- stringdist::stringdist(.x,index_list)
      # Remove those above barcode_mismatch threshold
      index_list <- index_list[index_dist <= barcode_mismatch]
      index_dist <- index_dist[index_dist <= barcode_mismatch]
      return(index_list[which.min(index_dist)])
    })) %>%
    tidyr::unnest(c(index, index2)) %>%
    group_by(index, index2, Sample_Name) %>%
    summarise(Freq = sum(Freq))

  gg.switch <- switch_plot_dat %>%
    dplyr::group_by(Sample_Name, index, index2) %>%
    summarise(Freq = sum(Freq))%>%
    dplyr::mutate(index = factor(index, levels=applied_indices$index), index2=factor(index2, levels=rev(applied_indices$index2)))  %>%
    ggplot2::ggplot(aes(x = index, y = index2), stat="identity") +
    geom_tile(aes(fill = Freq),alpha=0.8)  + 
    scale_fill_viridis_c(name="log10 Reads", begin=0.1, trans="log10")+
    theme(axis.text.x = element_text(angle=90, hjust=1), 
          plot.title=element_text(hjust = 0.5),
          plot.subtitle =element_text(hjust = 0.5)
    ) +
    labs(title= fcid, subtitle = paste0(
      "Total Reads: ", sum(indices$Freq, na.rm=TRUE),
      ", Switch rate: ", sprintf("%1.4f%%", res$switch_rate*100),
      ", Contam rate: ", sprintf("%1.6f%%", res$contam_rate*100),
      ", Other reads: ", other_reads)) 
  pdf(file=paste(fcid,"_index_switching.pdf"), width = 11, height = 8 , paper="a4r")
    plot(gg.switch)
    try(dev.off(), silent=TRUE)
  return(res)
}


# Plots -------------------------------------------------------------------
plot_read_quals <- function(sample_id, input_dir, truncLen = NULL, quiet=FALSE, n = 10000){
  input_dir <- normalizePath(input_dir)
  # Seq dir might need ot be changed to trimmed or other sub folders
  fastqFs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R1_*")))
  fastqRs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R2_*")))
  
  if(length(fastqFs) == 0 ){
    message(paste0("Sample ", sample_id, " Has no reads"))
    return(NULL)
  }
  if(!file.size(fastqFs) > 28) {
    message(paste0("Sample ", sample_id, " Has no reads"))
    return(NULL)
  }
  
  Fquals <- get_qual_stats(fastqFs, n=n)
  Rquals <- get_qual_stats(fastqRs, n=n)
  
  #Plot qualities
  gg.Fqual <- Fquals %>% 
    dplyr::select(Cycle, reads, starts_with("Q")) %>% 
    tidyr::pivot_longer(cols = starts_with("Q")) %>% 
    ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
    geom_line(data = Fquals, aes(y = QMean), color = "#66C2A5") + 
    geom_line(data = Fquals, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    geom_line(data = Fquals, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
    geom_line(data = Fquals, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
    ggtitle(paste0(sample_id, " Forward Reads")) +
    scale_x_continuous(breaks=seq(0,300,25))
    
  gg.Rqual <- Rquals %>% 
    dplyr::select(Cycle, reads, starts_with("Q")) %>% 
    tidyr::pivot_longer(cols = starts_with("Q")) %>% 
    ggplot2::ggplot(aes(x = Cycle, y = value, colour = name)) + 
    geom_line(data = Rquals, aes(y = QMean), color = "#66C2A5") + 
    geom_line(data = Rquals, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    geom_line(data = Rquals, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
    geom_line(data = Rquals, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    labs(x = "Reads position", y = "Quality Score") + ylim(c(0, NA))+
    ggtitle(paste0(sample_id, " Reverse Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")
  
  
  gg.Fee <- Fquals %>% 
    dplyr::select(Cycle, starts_with("EE")) %>% 
    tidyr::pivot_longer(cols = starts_with("EE")) %>% 
    dplyr::group_by(name) %>%  # Remove EE and replace with percentage at end? - i.e lower 10%
    dplyr::mutate(cumsumEE = cumsum(value)) %>%
    ggplot2::ggplot(aes(x = Cycle, y = log10(cumsumEE), colour = name)) + 
    geom_point(size = 1) + 
    geom_hline(yintercept = log10(1), color = "red") + 
    geom_hline(yintercept = log10(2), color = "red") + 
    geom_hline(yintercept = log10(3), color = "red") + 
    geom_hline(yintercept = log10(5), color = "red") + 
    geom_hline(yintercept = log10(7), color = "red") + 
    geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
    geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") + 
    labs(x = "Reads position", y = "Log10 Cumulative expected errors",
         colour = "Read quantiles") + 
    ggtitle(paste0(sample_id, " Forward Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")

  
  gg.Ree <- Rquals %>% 
    dplyr::select(Cycle, starts_with("EE")) %>% 
    tidyr::pivot_longer(cols = starts_with("EE")) %>% 
    dplyr::group_by(name) %>% 
    dplyr::mutate(cumsumEE = cumsum(value)) %>%
    ggplot2::ggplot(aes(x = Cycle, y = log10(cumsumEE), colour = name)) + 
    geom_point(size = 1) + 
    geom_hline(yintercept = log10(1), color = "red") + 
    geom_hline(yintercept = log10(2), color = "red") + 
    geom_hline(yintercept = log10(3), color = "red") + 
    geom_hline(yintercept = log10(5), color = "red") + 
    geom_hline(yintercept = log10(7), color = "red") + 
    geom_text(label = "MaxEE=1", aes(x = 0, y = log10(1), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=2", aes(x = 0, y = log10(2), hjust = 0, vjust = 0), color = "red") +
    geom_text(label = "MaxEE=3", aes(x = 0, y = log10(3), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=5", aes(x = 0, y = log10(5), hjust = 0, vjust = 0), color = "red") + 
    geom_text(label = "MaxEE=7", aes(x = 0, y = log10(7), hjust = 0, vjust = 0), color = "red") + 
    labs(x = "Reads position", y = "Log10 Cumulative expected errors") +
    ggtitle(paste0(sample_id, " Reverse Reads")) +
    scale_x_continuous(breaks=seq(0,300,25)) +
    theme(legend.position = "bottom")
  
  if(!is.null(truncLen)){
    gg.Fqual <- gg.Fqual +
      geom_vline(aes(xintercept=truncLen[1]), colour="blue") +
      annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    gg.Fee <-gg.Fee +
      geom_vline(aes(xintercept=truncLen[1]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[1]), colour="blue")
    gg.Rqual <- gg.Rqual +
      geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =2, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
    gg.Ree <- gg.Ree + 
      geom_vline(aes(xintercept=truncLen[2]), colour="blue")+
      annotate("text", x = truncLen[1]-10, y =-3, label = paste0("Suggested truncLen = ", truncLen[2]), colour="blue")
  }
  
  Qualplots <- (gg.Fqual + gg.Rqual) / (gg.Fee + gg.Ree)
  return(Qualplots)
}


# Primer trimming ---------------------------------------------------------
step_primer_trim <- function(sample_id, input_dir, output_dir, qc_dir, for_primer_seq, rev_primer_seq, pcr_primers,
                             max_mismatch = 0, n = 1e6, qualityType = "Auto", check_paired = TRUE, compress =TRUE, quiet=FALSE){
  input_dir <- normalizePath(input_dir)
  output_dir <- normalizePath(output_dir)
  qc_dir <- normalizePath(qc_dir)
  # Check inputs
  if(!is.character(for_primer_seq) | !length(for_primer_seq)==1){
    stop("for_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  if(!is.character(rev_primer_seq) | !length(rev_primer_seq)==1){
    stop("rev_primer_seq must be a character vector of length 2, containing forward and reverse primers")
  }
  # Check that required files exist
  if(!dir.exists(input_dir)) {
    stop("input_dir doesnt exist, check that the correct path was provided")
  }
  
  # Create output directory if it doesnt exist
  if(!dir.exists(output_dir)) {dir.create(output_dir)}
  
  # Create qc_dir if it doesnt exist
  if(!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}
  
  fastqFs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R1_*")))
  fastqRs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R2_*")))
  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))
  
  #Check if there were more than 1 primer per sample
  multi_primer <- any(stringr::str_detect(for_primer_seq, ";"), stringr::str_detect(rev_primer_seq, ";"), na.rm = TRUE)
  
  if (multi_primer) {
    for_primer_seq <- unlist(stringr::str_split(for_primer_seq, ";"))
    rev_primer_seq <- unlist(stringr::str_split(rev_primer_seq, ";"))
    
    primer_names <- unlist(stringr::str_split(pcr_primers, ";"))
    
    # Define outfiles
    fwd_out <- purrr::map(primer_names, ~{
      new_sampleid <- paste0(sample_id, "_",.x)
      return(normalizePath(paste0(output_dir,"/", stringr::str_replace(basename(fastqFs), sample_id, new_sampleid))))
    }) %>%
      unlist()
    rev_out <- purrr::map(primer_names, ~{
      new_sampleid <- paste0(sample_id, "_",.x)
      return(normalizePath(paste0(output_dir,"/", stringr::str_replace(basename(fastqRs), sample_id, new_sampleid))))
    }) %>%
      unlist()
    
  } else if (!multi_primer) {
    # CHange this to rename samples if primers not present
    new_sampleid <- paste0(sample_id, "_",pcr_primers)
    fwd_out <- normalizePath(paste0(output_dir,"/", stringr::str_replace(basename(fastqFs), sample_id, new_sampleid)))
    rev_out <- normalizePath(paste0(output_dir,"/", stringr::str_replace(basename(fastqRs), sample_id, new_sampleid)))
  } 
  
  # Demultiplex reads and trim primers
  res <- trim_primers(fwd = fastqFs,
                      rev = fastqRs,
                      fwd_out = fwd_out,
                      rev_out = rev_out,
                      for_primer_seq = stringr::str_replace_all(for_primer_seq, "I", "N"),
                      rev_primer_seq = stringr::str_replace_all(rev_primer_seq, "I", "N"),
                      n = n, qualityType = qualityType, check_paired = check_paired,
                      compress =compress, quiet=quiet
  ) %>%
    dplyr::mutate(fwd_out = fwd_out,
                  rev_out = rev_out)
  return(res)
}

# Primer trimming function
trim_primers <- function(fwd, rev, fwd_out, rev_out, for_primer_seq, rev_primer_seq, 
                         n = 1e6, qualityType = "Auto", check_paired = TRUE, id.field = NULL, 
                         max_mismatch=0, id.sep = "\\s", compress =TRUE, quiet=FALSE){
  
  ## iterating through forward and reverse files using fastq streaming
  fF <- ShortRead::FastqStreamer(file.path(fwd), n = n)
  on.exit(close(fF))
  fR <- ShortRead::FastqStreamer(file.path(rev), n = n)
  on.exit(close(fR), add=TRUE)
  
  #Check if there were more than 1 primer per sample
  if(any(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq) >0)){
    multi_primer <- TRUE
  }
  
  # Check if number of F and R primers match the number of outfiles
  if (!all.equal(length(fwd_out), length(rev_out), length(for_primer_seq), length(rev_primer_seq))) {
    stop("fwd_out and rev_out must be the same length as for_primer_seq and rev_primer_seq")
  }
  
  #if(!C_isACGT(primer)) stop("Non-ACGT characters detected in primers")
  
  # Delete output files if they already exist
  c(fwd_out, rev_out) %>%
    purrr::walk(remove_if_exists, quiet=quiet)
  
  first=TRUE
  append <- vector("logical", length= length(fwd_out)) 
  remainderF <- ShortRead::ShortReadQ(); remainderR <- ShortRead::ShortReadQ()
  casava <- "Undetermined" #ID field format
  # Setup read tracking
  inseqs <- 0
  outseqs <- vector("numeric", length= length(fwd_out))
  
  while( TRUE ) {
    suppressWarnings(fqF <- ShortRead::yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- ShortRead::yield(fR, qualityType = qualityType))
    if(length(fqF) == 0 && length(fqR) == 0) { break } # Stop loop if theres no reads left to process
    
    inseqs <- inseqs + length(fqF)
    
    # Make sure that forward and reverse reads are correctly paired
    # Determine the sequence identifier field. Looks for a single 6-colon field (CASAVA 1.8+ id format)
    if(check_paired) {
      if(first) { 
        if(is.null(id.field)) {
          # or a single 4-colon field (earlier format). Fails if it doesn't find such a field.
          id1 <- as.character(ShortRead::id(fqF)[[1]])
          id.fields <- strsplit(id1, id.sep)[[1]]
          ncolon <- sapply(gregexpr(":", id.fields), length)
          ncoltab <- table(ncolon)
          if(max(ncolon) == 6 && ncoltab["6"] == 1) { # CASAVA 1.8+ format
            casava <- "Current"
            id.field <- which(ncolon == 6)
          } else if (max(ncolon) == 4 && ncoltab["4"] == 1) { # CASAVA <=1.7 format
            casava <- "Old"
            id.field <- which(ncolon == 4)
          } else { # Couldn't unambiguously find the seq id field
            stop("Couldn't automatically detect the sequence identifier field in the fastq id string.")
          }
        }
      } else { # !first
        # Prepend the unmatched sequences from the end of previous chunks
        # Need ShortRead::append or the method is not dispatched properly
        fqF <- ShortRead::append(remainderF, fqF)
        fqR <- ShortRead::append(remainderR, fqR)
      }
    } else { # !check_paired
      if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files: ", length(fqF), ", ", length(fqR), ".")
    }
    
    # Enforce id matching (ASSUMES SAME ORDERING IN F/R, BUT ALLOWS DIFFERENTIAL MEMBERSHIP)
    # Keep the tail of unmatched sequences (could match next chunk)
    if(check_paired) {
      idsF <- sapply(strsplit(as.character(ShortRead::id(fqF)), id.sep), `[`, id.field)
      idsR <- sapply(strsplit(as.character(ShortRead::id(fqR)), id.sep), `[`, id.field)
      if(casava == "Old") { # Drop the index number/pair identifier (i.e. 1=F, 2=R)
        idsF <- sapply(strsplit(idsF, "#"), `[`, 1)
      }
      lastF <- max(c(0,which(idsF %in% idsR)))
      lastR <- max(c(0,which(idsR %in% idsF)))
      if(lastF < length(fqF)) {
        remainderF <- fqF[(lastF+1):length(fqF)]
      } else {
        remainderF <- ShortRead::ShortReadQ() 
      }
      if(lastR < length(fqR)) {
        remainderR <- fqR[(lastR+1):length(fqR)]
      } else {
        remainderR <- ShortRead::ShortReadQ() 
      }
      fqF <- fqF[idsF %in% idsR]
      fqR <- fqR[idsR %in% idsF]
    }
    
    # Demultiplex and trim each primer
    # Only keep reads where primer is detected
    
    for (p in 1:length(for_primer_seq)){
      barlenF <- nchar(for_primer_seq[p])
      barlenR <- nchar(rev_primer_seq[p])
      
      keepF <- Biostrings::neditStartingAt(
        pattern=Biostrings::DNAString(for_primer_seq[p]),
        subject= IRanges::narrow(sread(fqF), 1, barlenF),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE ) <= max_mismatch
      
      keepR <- Biostrings::neditStartingAt(
        pattern= Biostrings::DNAString(rev_primer_seq[p]),
        subject= IRanges::narrow(sread(fqR), 1, barlenR),
        starting.at=1,
        with.indels=FALSE,
        fixed=FALSE ) <= max_mismatch
      
      # Only keep reads where forward primer is present in F, and reverse in R
      keep <- keepF & keepF
      
      fqF_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqF[keep]), 
                                                           quality=Biostrings::quality(Biostrings::quality(fqF[keep])),
                                                           id=ShortRead::id(fqF[keep])))
      
      fqR_primer <- suppressWarnings(ShortRead::ShortReadQ(sread=ShortRead::sread(fqR[keep]), 
                                                           quality=Biostrings::quality(Biostrings::quality(fqR[keep])),
                                                           id=ShortRead::id(fqR[keep])))
      
      # Trim primers from left side
      startF <- max(1, barlenF + 1, na.rm=TRUE)
      startR <- max(1, barlenR + 1, na.rm=TRUE)
      
      # Make sure reads are longer than the primer sequence
      keep <- (width(fqF_primer) >= startF & width(fqR_primer) >= startR)
      fqF_primer <- fqF_primer[keep]
      fqF_primer <- narrow(fqF_primer, start = startF, end = NA)
      fqR_primer <- fqR_primer[keep]
      fqR_primer <- narrow(fqR_primer, start = startR, end = NA)
      
      outseqs[p] <- outseqs[p] + length(fqF_primer)    
      
      if(!append[p]) {
        ShortRead::writeFastq(fqF_primer, fwd_out[p], "w", compress = compress)
        ShortRead::writeFastq(fqR_primer, rev_out[p], "w", compress = compress)
        append[p]=TRUE
        first=FALSE
      } else {
        ShortRead::writeFastq(fqF_primer, fwd_out[p], "a", compress = compress)
        ShortRead::writeFastq(fqR_primer, rev_out[p], "a", compress = compress)
      }
    }
  }
  
  if(!quiet) {
    outperc <- purrr::map(outseqs, ~{
      round(.x * 100 / inseqs, 1)
    }) %>%
      unlist()
    outperc <- paste(" (", outperc, "%),", sep="")
    message("Read in ", inseqs, " paired-sequences, output ", paste(" ", outseqs, ",", sep=""), " ", outperc, " primer-trimmed paired-sequences.", sep="")
  }
  
  if(sum(outseqs)==0) {
    message(paste("No reads remaining for:", fwd, "and", rev))
    file.remove(fwd_out)
    file.remove(rev_out)
  }
  out <- data.frame(for_primer_seq = for_primer_seq,
                    rev_primer_seq = rev_primer_seq,
                    trimmed_input = inseqs,
                    trimmed_output = outseqs)
  return(out)
}


# Read filtering ----------------------------------------------------------

step_filter_reads <- function(sample_id, input_dir, output_dir, min_length = 20, max_length = Inf,
                              max_ee = 1, trunc_length = 150, trim_left = 0, trim_right = 0,
                              quiet=FALSE, ...){
  input_dir <- normalizePath(input_dir)
  output_dir <- normalizePath(output_dir)
  
  # Check that required files exist
  if(!dir.exists(input_dir)) {
    stop("input_dir doesnt exist, check that the correct path was provided")
  }
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {dir.create(output_dir)}
  
  fastqFs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R1_*")))
  fastqRs <- normalizePath(fs::dir_ls(path = input_dir, glob = paste0("*",sample_id, "_S*_R2_*")))

  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))

  # Handle NA inputs
  if(is.na(min_length)){min_length <- 20}
  if(is.na(max_length)){max_length <- Inf}
  if(is.na(max_ee)){max_ee <- Inf}
  if(is.na(trunc_length)){trunc_length <- 0}
  if(is.na(trim_left)){trim_left <- 0}
  if(is.na(trim_right)){trim_right <- 0}
  
  # Run read filter
  res <- dada2::filterAndTrim(
    fwd = fastqFs, filt = file.path(output_dir, basename(fastqFs)),
    rev = fastqRs, filt.rev = file.path(output_dir, basename(fastqRs)),
    minLen = min_length, maxLen = max_length, maxEE = max_ee, truncLen = trunc_length,
    trimLeft = trim_left, trimRight = trim_right, rm.phix = TRUE, 
    multithread = FALSE, compress = TRUE, verbose = !quiet) %>% 
    as_tibble() %>%
    dplyr::rename(filter_input = reads.in,
                  filter_output = reads.out) %>%
    dplyr::mutate(fwd_out = file.path(output_dir, basename(fastqFs)),
           rev_out = file.path(output_dir, basename(fastqRs)))

  if(res$filter_output > 0){
    filtered_summary <- res %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::mutate(reads_remaining = signif(((filter_output / filter_input) * 100), 2)) %>%
      dplyr::filter(!is.na(reads_remaining))
    
    if (filtered_summary$reads_remaining < 10) {
      message(paste0("WARNING: Less than 10% reads remaining for ",
                     filtered_summary$sample), "Check filtering parameters ")
    }
  }
  return(res)
}

#  DADA2 ------------------------------------------------------------------
step_errormodel <- function(fcid, input_dir, pcr_primers, output, qc_dir, read="F", nbases=1e+08, 
                            randomize=FALSE, multithread=FALSE, write_all = FALSE, quiet=FALSE){
  
  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  if(read == "F"){
    filts <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
    message(paste0("Modelling forward read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
  } else if (read == "R"){
    filts <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
    message(paste0("Modelling reverse read error rates for primers: ", pcr_primers, " and flowcell: ", fcid))
  } else {
    stop ("read must be F or R!")
  }
  # Subset fastqs to just the relevent pcr primers
  filts <- filts[str_detect(filts,paste0(pcr_primers, "(-|_|$)"))]
  message(paste0(length(filts), " fastq files to process for primers: ", pcr_primers, " and flowcell: ", fcid))
  
  # Learn error rates from a subset of the samples and reads (rather than running self-consist with full dataset)
  err <-  dada2::learnErrors(filts, multithread = multithread, nbases = nbases,
                              randomize = randomize, qualityType = "FastqQuality", verbose=TRUE)
  #write out errors for diagnostics
  if(write_all){
    write_csv(as.data.frame(err$trans), paste0(qc_dir, "/", fcid, "_err",read,"_observed_transitions.csv"))
    write_csv(as.data.frame(err$err_out), paste0(qc_dir, "/", fcid, "_err",read,"_inferred_errors.csv"))
  }
  
  ##output error plots to see how well the algorithm modelled the errors in the different runs
  p1 <- dada2::plotErrors(err, nominalQ = TRUE) + ggtitle(paste0(pcr_primers, " ", fcid, " Forward Reads"))
  pdf(paste0(qc_dir,"/",fcid,"_", pcr_primers, "_", read,"_errormodel.pdf"), width = 11, height = 8 , paper="a4r")
    plot(p1)
  try(dev.off(), silent=TRUE)
  
  saveRDS(err, output)
}

step_dada2 <- function(fcid, input_dir, pcr_primers, output, qc_dir, error_model, pool="pseudo",
                       quiet=FALSE,  multithread=FALSE){
  
  errF <- error_model[[1]]
  errR <- error_model[[2]]
  
  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  filtFs <- list.files(input_dir, pattern="R1_001.*", full.names = TRUE)
  filtRs <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
  
  # Subset fastqs to just the relevent pcr primers
  filtFs <- filtFs[str_detect(filtFs,paste0(pcr_primers, "(-|_|$)"))]
  filtRs <- filtRs[str_detect(filtRs,paste0(pcr_primers, "(-|_|$)"))]
  if(length(filtFs) != length(filtRs)) stop(paste0("Forward and reverse files for ",fcid," do not match."))
  message(paste0(length(filtFs), " fastq files to process for primers: ", pcr_primers, " and flowcell: ", fcid))
  
  #Denoise reads
  message(paste0("Denoising forward reads for primers: ", pcr_primers, " and flowcell: ", fcid))
  dadaFs <- dada2::dada(filtFs, err = errF, multithread = multithread, pool = pool, verbose = TRUE)
  message(paste0("Denoising reverse reads for primers: ", pcr_primers, " and flowcell: ", fcid))
  dadaRs <- dada2::dada(filtRs, err = errR, multithread = multithread, pool = pool, verbose = TRUE)
  
  dada <- list(dadaFs, dadaRs)
  saveRDS(dada, output)
}


step_dada2_single <- function(fcid, sample_id, input_dir, pcr_primers, output, qc_dir, error_model,
                              priors_fwd = NA, priors_rev = NA, quiet=FALSE,  multithread=FALSE){
  
  errF <- error_model[[1]]
  errR <- error_model[[2]]
  
  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  filtFs <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
  filtRs <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
  
  # Subset fastqs to just the relevent sample_id
  filtFs <- filtFs[str_detect(filtFs, sample_id)]
  filtRs <- filtRs[str_detect(filtRs, sample_id)]
  if(length(filtFs) != length(filtRs)) stop(paste0("Forward and reverse files for ",sample_id," do not match."))

  # Handle priors
  if(length(priors_fwd) == 1 && is.na(priors_fwd)){
    priors_fwd <- character(0)
  }
  if(length(priors_rev) == 1 && is.na(priors_rev)){
    priors_rev <- character(0)
  }
  #Denoise reads
  message(paste0("Denoising forward and reverse reads for sample: ", sample_id))
  if(length(priors_fwd) > 1 | length(priors_rev) > 1){
    message("High Sensitivity mode set: Using prior sequences to refine ASV inference")
  }
  dadaFs <- dada2::dada(filtFs, err = errF, multithread = multithread, priors = priors_fwd, selfConsist = FALSE, pool = FALSE, verbose = TRUE)
  dadaRs <- dada2::dada(filtRs, err = errR, multithread = multithread, priors = priors_rev, selfConsist = FALSE, pool = FALSE, verbose = TRUE)
  
  dada <- list(dadaFs, dadaRs)
  saveRDS(dada, output)
}

step_dada2_single2 <- function(fcid, sample_id, input_dir, pcr_primers, output, qc_dir, error_model, read="F",
                               priors = NA, quiet=FALSE,  multithread=FALSE){

  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  
  if(read == "F"){
    filts <- list.files(input_dir, pattern= "*R1_001.*", full.names = TRUE)
    message(paste0("Denoising forward reads for sample: ", sample_id))
  } else if (read == "R"){
    filts <- list.files(input_dir, pattern="R2_001.*", full.names = TRUE)
    message(paste0("Denoising reverse reads for sample: ", sample_id))
  } else {
    stop ("read must be F or R!")
  }
  # Subset fastqs to just the relevent sample_id
  filts <- filts[str_detect(filts, sample_id)]
  
  # Handle priors
  if(length(priors) == 1 && is.na(priors)){
    priors <- character(0)
  }
  #Denoise reads
  if(length(priors) > 1 ){
    message("High Sensitivity mode set: Using prior sequences to refine ASV inference")
  }
  dada <- dada2::dada(filts, err = error_model, multithread = multithread, priors = priors, selfConsist = FALSE, pool = FALSE, verbose = TRUE)
  saveRDS(dada, output)
}


step_mergereads <- function(fcid, input_dir, pcr_primers, output, qc_dir, dada,
                            write_all = FALSE, quiet=FALSE,  multithread=FALSE, concat_unmerged=FALSE){
  
  # read in denoised files and subset to just the relevent pcr primers
  dadaFs <- dada[[1]]
  dadaRs <- dada[[2]]
  dadaFs <- dadaFs[str_detect(names(dadaFs),paste0(pcr_primers, "(-|_|$)"))]
  dadaRs <- dadaRs[str_detect(names(dadaRs),paste0(pcr_primers, "(-|_|$)"))]
  
  input_dir <- normalizePath(input_dir)
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  
  # Read in fastqs for just those samples in the dadas
  filtFs <- map_chr(names(dadaFs), ~{
    list.files(input_dir, pattern=paste0(.x, ".*R1_001.*"), full.names = TRUE)
  })
  filtRs <- map_chr(names(dadaRs), ~{
    list.files(input_dir, pattern=paste0(.x, ".*R2_001.*"), full.names = TRUE)
  })

  if(!all.equal(length(filtFs),length(filtRs), length(dadaFs), length(dadaFs))){
    stop(paste0("Number of input files dont match! (filtered F:",
                length(filtFs), ", filtered R: ", length(filtFs), 
                ", denoised F:", length(dadaFs), ", denoised R:", length(dadaRs),")"))
    
  }

  # Merge reads
  message(paste0("Merging forward and reverse reads for ", pcr_primers, " and flowcell: ", fcid))
  if(write_all | concat_unmerged){
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 12, trimOverhang = TRUE, returnRejects=TRUE) 
  } else {
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 12, trimOverhang = TRUE) 
  }
  
  # Write out unmerged reads
  if(write_all){
    unmergedFs <- lapply(1:length(mergers), function(x) {
      if(any(!mergers[[x]]$accept)){
        tmp <- dadaFs[[x]]$sequence[mergers[[x]]$forward[!mergers[[x]]$accept]]
        names(tmp) <- paste0(x, "_", which(!mergers[[x]]$accept), "/1")
        return(tmp)
      } else {
        return(NULL)
      }
    }
    )
    unmergedRs <- lapply(1:length(mergers), function(x) {
      if(any(!mergers[[x]]$accept)){
        tmp <- dadaRs[[x]]$sequence[mergers[[x]]$reverse[!mergers[[x]]$accept]]
        names(tmp) <- paste0(x, "_", which(!mergers[[x]]$accept), "/2")
        return(tmp)
      } else {
        return(NULL)
      }
    }
    )
    # Write out fasta of concatenated reads
    unmerged_fasta <- paste0(unlist(unmergedFs), "NNNNNNNNNN", dada2::rc(unlist(unmergedRs)))
    names(unmerged_fasta) <- unlist(sapply(unmergedFs, names)) %>% stringr::str_remove("/.*$")
    taxreturn::write_fasta(taxreturn::char2DNAbin(unmerged_fasta), file = paste0(qc_dir, "/", fcid, "_unmerged_asvs.fa"))
  }
  
  # concatenate unmerged reads
  #Modified from https://github.com/benjjneb/dada2/issues/565
  if(concat_unmerged){
    message("concat_unmerged is set to TRUE - Concatenating unmerged forward and reverse reads")
    mergers_rescued <- mergers
    for(i in 1:length(mergers)) {
      if(any(!mergers[[i]]$accept)){
        # Get index of unmerged reads in table
        unmerged_index <- which(!mergers[[i]]$accept)
        # Get the forward and reverse reads for those reads
        unmerged_fwd <- dadaFs[[i]]$sequence[mergers[[i]]$forward[unmerged_index]]
        unmerged_rev <- dadaRs[[i]]$sequence[mergers[[i]]$reverse[unmerged_index]]
        
        unmerged_concat <- paste0(unmerged_fwd, "NNNNNNNNNN", rc(unmerged_rev))
        
        mergers_rescued[[i]]$sequence[unmerged_index] <- unmerged_concat
        mergers_rescued[[i]]$nmatch[unmerged_index] <- 0
        mergers_rescued[[i]]$nmismatch[unmerged_index] <- 0
        mergers_rescued[[i]]$nindel[unmerged_index] <- 0
        mergers_rescued[[i]]$prefer[unmerged_index] <- NA
        mergers_rescued[[i]]$accept[unmerged_index] <- TRUE
      } 
      
    }
    mergers <- mergers_rescued
  }
  
  if(write_all){
    dplyr::bind_rows(mergers, .id="Sample") %>%
      dplyr::mutate(Sample = stringr::str_remoce(Sample, pattern="_S[0-9]+_R[1-2]_.*$")) %>%
      write_csv(paste0(qc_dir, "/", fcid, "_mergers.csv"))
  }
  
  #Construct sequence table
  sample_names <- basename(filtFs) %>% stringr::str_remove("_S[0-9]+_R[1-2]_.*$")
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, output)
  
  # Track reads
  getN <- function(x) sum(getUniques(x))
  res <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN)) %>%
    magrittr::set_colnames(c("dadaFs", "dadaRs", "merged")) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    dplyr::mutate(sample_id = stringr::str_remove(basename(sample_id), pattern="_S[0-9]+_R[1-2]_.*$")) %>%
    as_tibble()
  return(res)
}


# ASV filtering -----------------------------------------------------------

# Group by target loci and apply
step_filter_asvs <- function(seqtab, pcr_primers, output, qc_dir, min_length = NULL, max_length = NULL, 
                             check_frame=FALSE, genetic_code="SGC4", phmm=NULL,
                             primers=NULL, multithread=FALSE, quiet=FALSE){
  
  # Handle NA inputs
  if(is.na(min_length)){min_length <- NULL}
  if(is.na(max_length)){max_length <- NULL}
  if(is.na(phmm)){phmm <- NULL}
  if(any(is.na(primers))){primers <- NULL}
  
  # Normalise paths
  output <- normalizePath(output)
  qc_dir <- normalizePath(qc_dir)
  
  if(is.matrix(seqtab) | is.data.frame(seqtab)){
    if(!quiet){message("Input is a matrix or data frame")}
  } else if (is.character(seqtab) & stringr::str_detect(seqtab, ".rds")){
    seqtab <- readRDS(seqtab)
  } else {
    stop("seqtab must be a matrix/data frame or .rds file")
  }
  reads_starting <- rowSums(seqtab)
  
  # Load in profile hidden markov model if provided
  if(is.character(phmm) && stringr::str_detect(phmm, ".rds")){
    phmm_model <- readRDS(phmm)
  } else if (is(phmm, "PHMM")){
    phmm_model <- phmm
  } else {
    phmm_model <- NULL
  }
  
  # subset PHMM if primers were provided
  if (is(phmm_model, "PHMM") && !is.null(primers)){
    # Check that one of the two primers can bind
    Fbind <- get_binding_position(primers[1], model = phmm_model, tryRC = TRUE, min_score = 10)
    Rbind <- get_binding_position(primers[2], model = phmm_model, tryRC = TRUE, min_score = 10)
    if(!is.na(Fbind$start) & !is.na(Rbind$start)){
      phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else  if(!is.na(Fbind$start) & is.na(Rbind$start)){
      # Reverse primer not found - Try with subsets
      for(r in seq(1, nchar(primers[2])-10, 1)){ #Minimum length of 10 as this has to match minscore
        Rbind <- get_binding_position(str_remove(primers[2], paste0("^.{1,",r,"}")), model = phmm_model, tryRC = TRUE, min_score = 10)
        if (!is.na(Rbind$start)) {
          primers[2] <- str_remove(primers[2], paste0("^.{1,",r,"}"))
          break
        }
      }
      phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    } else  if(is.na(Fbind$start) & !is.na(Rbind$start)){
      # Forward primer not found - Try with subsets
      for(r in seq(1, nchar(primers[1])-10, 1)){ #Minimum length of 10 as this has to match minscore
        Rbind <- get_binding_position(str_remove(primers[1], paste0("^.{1,",r,"}")), model = phmm_model, tryRC = TRUE, min_score = 10)
        if (!is.na(Rbind$start)) {
          primers[1] <- str_remove(primers[1], paste0("^.{1,",r,"}"))
          break
        }
      }
      phmm_model <- taxreturn::subset_model(phmm_model, primers = primers)
    }
  }
  # Remove chimeras  
  seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=!quiet)
  seqs_rem <- length(colnames(seqtab_nochim))/length(colnames(seqtab))
  abund_rem <- sum(seqtab_nochim)/sum(seqtab)
  message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% initial abundance remaining after chimera removal"))
  reads_chimerafilt <- rowSums(seqtab_nochim)
  
  # cut to expected size
  if(any(!is.null(c(min_length, max_length)), na.rm = TRUE) & any(reads_chimerafilt > 0)){
    if(!is.null(min_length) & !is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) %in% min_length:max_length]
    } else if(is.null(min_length) & !is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) < max_length]
    } else if(!is.null(min_length) & is.null(max_length)){
      seqtab_cut <- seqtab_nochim[,nchar(colnames(seqtab_nochim)) > min_length]
    }
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_cut))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_cut)/sum(seqtab)
      message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% initial abundance remaining after length filtering"))
    }
    reads_lengthfilt <- rowSums(seqtab_cut)
  } else {
    seqtab_cut <- seqtab_nochim
    reads_lengthfilt <- rep(0,nrow(seqtab)) 
    names(reads_lengthfilt) <- rownames(seqtab)
  }
  
  # Align against phmm
  if (is(phmm_model, "PHMM") & any(reads_lengthfilt > 0)){
    seqs <- DNAStringSet(colnames(seqtab_cut))
    names(seqs) <- colnames(seqtab_cut)
    phmm_filt <- taxreturn::map_to_model(
      seqs, model = phmm_model, min_score = 100, min_length = 100,
      shave = FALSE, check_frame = check_frame, kmer_threshold = 0.5, k=5, extra = "fill")
    seqtab_phmm <- seqtab_cut[,colnames(seqtab_cut) %in% names(phmm_filt)]
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_phmm))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_phmm)/sum(seqtab)
      message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% of initial abundance remaining after PHMM filtering"))
    }
    reads_phmmfilt <- rowSums(seqtab_phmm)
  } else {
    seqtab_phmm <- seqtab_cut
    reads_phmmfilt <- rep(0,nrow(seqtab)) 
    names(reads_phmmfilt) <- rownames(seqtab)
  }
  
  #Filter sequences containing stop codons
  if(check_frame & any(reads_phmmfilt > 0)){
    seqs <- DNAStringSet(colnames(seqtab_phmm))
    names(seqs) <- colnames(seqtab_phmm)
    codon_filt <- taxreturn::codon_filter(seqs, genetic_code = genetic_code) 
    seqtab_final <- seqtab_phmm[,colnames(seqtab_phmm) %in% names(codon_filt)]
    if(!quiet){
      seqs_rem <- length(colnames(seqtab_final))/length(colnames(seqtab))
      abund_rem <- sum(seqtab_final)/sum(seqtab)
      message(paste(round(seqs_rem*100,2), "% of initial sequences and ",round(abund_rem*100,2), "% of initial abundance remaining after checking reading frame"))
    }
    reads_framefilt <- rowSums(seqtab_final)
  } else {
    seqtab_final <- seqtab_phmm
    reads_framefilt <- rep(0,nrow(seqtab))
    names(reads_framefilt) <- rownames(seqtab)
  }
  reads_final <- rowSums(seqtab_final)
  
  saveRDS(seqtab_final, output)
        
  # Output a cleanup summary
  cleanup <- seqtab %>%
    as.data.frame() %>%
    pivot_longer( everything(),
                  names_to = "OTU",
                  values_to = "Abundance") %>%
    dplyr::group_by(OTU) %>%
    summarise(Abundance = sum(Abundance)) %>%
    dplyr::mutate(length  = nchar(OTU)) %>%
    dplyr::mutate(type = case_when(
      !OTU %in% getSequences(seqtab_nochim) ~ "Chimera",
      !OTU %in% getSequences(seqtab_cut) ~ "Incorrect size",
      !OTU %in% getSequences(seqtab_phmm) ~ "PHMM",
      !OTU %in% getSequences(seqtab_final) ~ "Stop codons",
      TRUE ~ "Retained"
    )) %>%
    dplyr::mutate(concat = str_detect(OTU, "NNNNNNNNNN")) %>%
    dplyr::mutate(type = case_when(
      concat ~ paste0("Unmerged-", type),
      TRUE ~ type
    )) %>%
    dplyr::mutate(pcr_primers = pcr_primers) %>%
    dplyr::select(-concat)
  
  cols <- c(`Chimera` = "#9e0142",
            `Unmerged-Chimera` = "#d53e4f",
            `Incorrect size` = "#f46d43",
            `Unmerged-Incorrect size` = "#fdae61",
            `PHMM` = "#fee08b",
            `Unmerged-PHMM` = "#e6f598",
            `Stop codons` = "#abdda4",
            `Unmerged-Stop codons` = "#66c2a5",
            `Retained` = "#3288bd",
            `Unmerged-Retained` = "#5e4fa2") 
  
  # Output length distribution plots
  gg.abundance <- ggplot2::ggplot(cleanup, aes(x=length, y=log10(Abundance), fill=type))+
    geom_bar(stat="identity") + 
    scale_x_continuous(limits=c(min(cleanup$length)-10, max(cleanup$length)+10))+
    theme_bw()+
    scale_fill_manual(values = cols)+
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
      title=pcr_primers,
      subtitle = "Abundance of sequences",
      x = "ASV length",
      y = "log10 ASV abundance",
      fill = "ASV type")
  
  gg.unique <- ggplot2::ggplot(cleanup, aes(x=length, fill=type))+
    geom_histogram(binwidth = 1) + 
    scale_x_continuous(limits=c(min(cleanup$length)-10, max(cleanup$length)+10))+
    theme_bw()+
    scale_fill_manual(values = cols)+
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
      title=pcr_primers,
      subtitle = "Number of unique sequences",
      x = "ASV length",
      y = "Number of unique sequences",
      fill = "ASV type")
  
  # Create combined plot
  out_plot <- gg.abundance / gg.unique
  
  # Create output
  res <- tibble(
    sample_id = rownames(seqtab) %>% stringr::str_remove(pattern="_S[0-9]+_R[1-2]_.*$"),
    reads_starting = reads_starting,
    reads_chimerafilt = reads_chimerafilt,
    reads_lengthfilt = reads_lengthfilt,
    reads_phmmfilt = reads_phmmfilt,
    reads_framefilt = reads_framefilt,
    reads_final = reads_final
  )
  return(list(filtered_seqtab = seqtab_final, 
              filtered_asvs = res,
              cleanup_summary = cleanup,
              plot = list(out_plot)))
}


# Taxonomic assignment ----------------------------------------------------

# Group by target loci and apply
step_idtaxa <- function(seqtab, qc_dir, database, threshold = 60, multithread=FALSE, quiet=FALSE,
                        ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                        return_ids=FALSE, remove_Ns=FALSE){
  # Load the relevent db
  trainingSet <- readRDS(normalizePath(database))
  
  # get the sequences from the seqtab
  seqs <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs
  # Stop if seqs are 0
  if(length(seqs) > 0){
    
    # Remove any 10bp N bases that were added by concatenating reads reads  
    if(remove_Ns){
      if(any(seqs %>% purrr::map_lgl(~{str_detect(as.character(.x), "NNNNNNNNNN")}))){
        seqs <- DNAStringSet(seqs %>% purrr::map_chr(~{str_replace(as.character(.x), "NNNNNNNNNN", "")}))
      }
    }
    
    # Classify 
    ids <- DECIPHER::IdTaxa(seqs, trainingSet, processors=1, threshold = threshold, verbose=!quiet, strand = "top") 
    
    # Get the filename of that db that we can use to name the output files
    db_name <- basename(database) %>% stringr::str_remove("\\..*$") %>% stringr::str_remove("_idtaxa")
    
    # Check that more than just root has been assigned
    if( any(sapply(ids, function(x){ length(x$taxon) }) > 2)){
      #Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
      tax <- ids %>%
        purrr::map_dfr(function(x){
          taxa <- paste0(x$taxon,"_", x$confidence)
          taxa[startsWith(taxa, "unclassified_")] <- NA
          data.frame(t(taxa)) %>%
            magrittr::set_colnames(ranks[1:ncol(.)])
        }) %>%
        mutate_all(stringr::str_replace,pattern="(?:.(?!_))+$", replacement="") %>%
        magrittr::set_rownames(getSequences(seqtab))
      # add empty ranks if none were assigned to lower ranks
      tax <- new_bind(tibble::tibble(!!!ranks, .rows = 0, .name_repair = ~ ranks), tax)
    } else {
      warning(paste0("No sequences assigned with IDTAXA to ", database, " have you used the correct database?"))
      tax <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") 
      tax[ranks[1]] <- ranks[1]
      tax[ranks[2:length(ranks)]] <- NA_character_
      tax <- tax %>%
        magrittr::set_rownames(getSequences(seqtab)) 
    }
  } else {
    warning(paste0("No sequences present in seqtab - IDTAXA skipped"))
    tax <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") 
    tax[ranks[1]] <- ranks[1]
    tax[ranks[2:length(ranks)]] <- NA_character_
    tax <- tax %>%
      magrittr::set_rownames(getSequences(seqtab)) 
  }
  # Check that output dimensions match input
  if(!all(rownames(tax) %in% colnames(seqtab))){
    stop("Number of ASVs classified does not match the number of input ASVs")
  }
  
  # Check that all ranks are present
  if(!all(colnames(tax) %in% ranks)){
    stop("Number of ranks does not match")
  }
  
  # Return objects
  if(return_ids){
    out <- list(tax = tax,
                ids = ids)
  } else {
    out <- tax
  }
  return(out)
}

step_blast_tophit <- function(seqtab, output=NULL, qc_dir, database, identity = 97,  coverage=95, evalue=1e06,
                              max_target_seqs=5, max_hsp=5, 
                              ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                              multithread=FALSE, quiet=FALSE){
  
  seqmap <- tibble::enframe(getSequences(seqtab), name = NULL, value="OTU") %>%
    mutate(name = paste0("SV", seq(length(getSequences(seqtab)))))
  
  seqs <- taxreturn::char2DNAbin(seqmap$OTU)
  names(seqs) <- seqmap$name
  
  # Get the filename of that db that we can use to name the output files
  db_name <- basename(database) %>% stringr::str_remove("_.*$")
  
  # Stop if seqs are 0
  if(length(seqs) > 0){
    blast_spp <- blast_assign_species(query=seqs,db=database, identity=97, coverage=95, evalue=1e06,
                                      max_target_seqs=5, max_hsp=5, ranks=ranks, delim=";") %>%
      dplyr::rename(blast_genus = Genus, blast_spp = Species) %>%
      dplyr::filter(!is.na(blast_spp)) 
    
    if(nrow(blast_spp) > 0){
      # Transform into taxtab format
      out <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") %>%
        dplyr::left_join(blast_spp %>%
                           dplyr::select(name = OTU, Genus = blast_genus, Species = blast_spp) %>%
                           left_join(seqmap) %>%
                           dplyr::select(-name), by="OTU")%>%
        column_to_rownames("OTU") %>%
        as.matrix()
    } else{
      warning(paste0("No Species assigned with BLAST to ", database, " have you used the correct database?"))
      out <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") %>%
        dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
        column_to_rownames("OTU") %>%
        as.matrix()
    }
  } else {
    warning(paste0("No sequences present in seqtab - BLAST skipped"))
    out <- tibble::enframe(getSequences(seqtab), name=NULL, value="OTU") %>%
      dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
      column_to_rownames("OTU") %>%
      as.matrix()
  }
  
  # Check that output dimensions match input
  if(!all(rownames(out) %in% colnames(seqtab))){
    stop("Number of ASVs classified does not match the number of input ASVs")
  }
  
  if(!is.null(output)){saveRDS(out, output)}
  return(out)
}

coalesce_tax <- function (x, y, suffix = c(".x", ".y"), prefer="left", join = dplyr::full_join, 
                          ...) {
  if(!"OTU" %in% colnames(x)){
    x <- x %>%
      as_tibble(rownames = "OTU")
  }
  if(!"OTU" %in% colnames(y)){
    y <- y %>%
      as_tibble(rownames = "OTU")
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
    column_to_rownames("OTU")
  return(out)
}


step_join_tax_blast <- function(tax, blast_spp, output=NULL, propagate_tax=FALSE,
                          ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")){
  if(!all(rownames(tax) %in% rownames(blast_spp))){
    stop("ASVs in tax and blast_spp do not match")
  }
  
  if(nrow(tax) > 0 & nrow(blast_spp) > 0){
    #Set BLAST ids where tax isnt at genus level to NA
    disagreements <- !blast_spp[,1] == tax[,7]
    disagreements[is.na(disagreements)] <- TRUE
    blast_spp[,1][disagreements] <- NA
    blast_spp[,2][disagreements] <- NA
    
    # Join the two taxonomies, prefering names from the tax
    tax_blast <- coalesce_tax(tax, blast_spp)
  
    if(propagate_tax){
      tax_blast <- tax_blast %>%
        seqateurs::na_to_unclassified() #Propagate high order ranks to unassigned ASVs
    }
  } else {
    warning("Either tax or blast_spp is empty")
    tax_blast <- as.matrix(tax)
  }
  if(!all(rownames(tax) %in% rownames(tax_blast))){
    stop("Number of ASVs output does not match the number of input ASVs")
  }
  # Write taxonomy table for that db to disk
  if(!is.null(output)){saveRDS(tax_blast, output)}
  return(tax_blast)
}


# Outputs -----------------------------------------------------------------

step_phangorn <- function(seqtab, output=NULL){

  seqs <- getSequences(seqtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  
  #Fit NJ tree
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit <- pml(treeNJ, data=phang.align)
  
  #Fit ML tree
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  
  # Write phytree to disk
  saveRDS(fitGTR, "output/rds/phytree.rds") 
  
  #Output newick tree
  if(!is.null(output)){write.tree(fitGTR$tree, file=output)}
  
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
      samdf <- read_csv(normalizePath(samdf))
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
    seqs <- DNAStringSet(colnames(seqtab))
    names(seqs) <- seqs
  }
  
  #Check if phy is a path - if so read in
  if(is(phylo, "character")){
    if (file.exists(phylo)){
      phy <- read.tree(normalizePath(phylo))
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
    filter(!duplicated(sample_id)) %>%
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
    taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)),"-",tax_table(ps)[,8])
  }
  
  if(nrow(seqtab) > nrow(phyloseq::sample_data(ps))){
    message("Warning: the following samples were not included in phyloseq object, check sample names match the sample metadata")
    message(rownames(seqtab)[!rownames(seqtab) %in% sample_names(ps)])
  }
  return(ps)
}

step_rareplot <- function(ps, min_reads=1000, plot_dir=NULL){
  #Create rarefaction curve
  rare <- phyloseq::otu_table(ps) %>%
    as("matrix") %>%
    rarecurve(step=max(sample_sums(ps))/100) %>%
    purrr::map(function(x){
      b <- as.data.frame(x)
      b <- data.frame(OTU = b[,1], count = rownames(b))
      b$count <- as.numeric(gsub("N", "",  b$count))
      return(b)
    }) %>%
    purrr::set_names(sample_names(ps)) %>%
    dplyr::bind_rows(.id="sample_id")
  
  gg.rare <- ggplot2::ggplot(data = rare)+
    geom_line(aes(x = count, y = OTU, group=sample_id), alpha=0.5)+
    geom_point(data = rare %>% 
                 dplyr::group_by(sample_id) %>% 
                 top_n(1, count),
               aes(x = count, y = OTU, colour=(count > min_reads))) +
    geom_label(data = rare %>% 
                 dplyr::group_by(sample_id) %>% 
                 top_n(1, count),
               aes(x = count, y = OTU,label=sample_id, colour=(count > min_reads)),
               hjust=-0.05)+
    scale_x_continuous(labels =  scales::scientific_format()) +
    geom_vline(xintercept=min_reads, linetype="dashed") +
    labs(colour = "Sample kept?") +
    xlab("Sequence reads") +
    ylab("Observed ASV's")
  
  #Write out figure
  pdf(file=normalizePath(paste0(plot_dir, "/rarefaction.pdf")), width = 11, height = 8 , paper="a4r")
  plot(gg.rare)
  try(dev.off(), silent=TRUE)
  
}

step_filter_phyloseq <- function(ps, kingdom = NA, phylum = NA, class = NA, 
                                 order = NA, family = NA, genus = NA, species = NA,
                                 min_sample_reads=1000, min_taxa_reads=NA, min_taxa_ra = NA, quiet=FALSE){
  
  #Taxonomic filtering
  taxtab <- phyloseq::tax_table(ps) %>%
    as("matrix") %>% 
    as.data.frame()
  
  # Check if any taxonomic filters are enabled
  ps0 <- ps
  if(any(!sapply(c(kingdom, phylum, class, order, family, genus, species), is.na))){
    
    # Filter kingdom
    if(!is.na(kingdom)){
      if (any(stringr::str_detect(taxtab$Kingdom, kingdom))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Kingdom",
            value = kingdom
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Kingdom", kingdom," - Check your target_kingdom parameter"))
      }
    }
    # Filter phylum
    if(!is.na(phylum)){
      if (any(stringr::str_detect(taxtab$Phylum, phylum))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Phylum",
            value = phylum
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Phylum", phylum," - Check your target_phylum parameter"))
      }
    }
    # Filter class
    if(!is.na(class)){
      if (any(stringr::str_detect(taxtab$Class, class))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Class",
            value = class
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Class", class," - Check your target_phylum parameter"))
      }
    }
    # Filter order
    if(!is.na(order)){
      if (any(stringr::str_detect(taxtab$Order, order))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Order",
            value = order
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Order", order," - Check your target_phylum parameter"))
      }
    }
    # Filter family
    if(!is.na(family)){
      if (any(stringr::str_detect(taxtab$Family, family))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Family",
            value = family
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Family", family," - Check your target_phylum parameter"))
      }
    }
    # Filter genus
    if(!is.na(family)){
      if (any(stringr::str_detect(taxtab$Genus, genus))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Genus",
            value = genus
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Genus", genus," - Check your target_phylum parameter"))
      }
    }
    # Filter Species
    if(!is.na(family)){
      if (any(stringr::str_detect(taxtab$Species, species))){
        ps0 <- ps0 %>%
          subset_taxa_new(
            rank = "Species",
            value = species
          ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
      } else{
        warning(paste0("No ASVs were assigned to the Species", species," - Check your target_phylum parameter"))
      }
    }
  } else {
    if (!quiet){message(paste0("No taxonomic filters set - skipping this filter"))}
  }
  
  # Remove any taxa under read count or relative abundance thresholds
  if(!is.na(min_taxa_reads) & is.na(min_taxa_ra)){
    ps1 <- phyloseq::transform_sample_counts(ps0, function(OTU, ab = min_taxa_reads){ ifelse(OTU <= ab,  0, OTU) })
  } else if(is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps1 <- phyloseq::transform_sample_counts(ps0, function(OTU, ab = min_taxa_ra ){
      ifelse((OTU / sum(OTU)) <= ab,  0, OTU) 
    })
  } else if (!is.na(min_taxa_reads) & !is.na(min_taxa_ra)){
    ps1 <- ps0 %>%
      phyloseq::transform_sample_counts(function(OTU, ab = min_taxa_reads){ ifelse(OTU <= ab,  0, OTU) }) %>%
      phyloseq::transform_sample_counts(function(OTU, ab = min_taxa_ra ){
        ifelse((OTU / sum(OTU)) <= ab,  0, OTU) 
      })
  } else {
    if (!quiet){message(paste0("No minimum abundance filters set - skipping this filter"))}
    ps1 <- ps0
  }
  
  #Remove all samples under the minimum read threshold 
  if(min_sample_reads > 0){
    ps2 <- ps1 %>%
      prune_samples(sample_sums(.)>=min_sample_reads, .) %>% 
      filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
  } else {
    if (!quiet){message(paste0("No minimum sample reads filter set - skipping this filter"))}
    ps2 <- ps1 %>%
      prune_samples(sample_sums(.)>=0,.) %>%
      filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
  }
  
  #Message how many were removed
  if(!quiet){message(nsamples(ps) - nsamples(ps2), " Samples and ", ntaxa(ps) - ntaxa(ps2), " ASVs dropped")}
  return(ps2)
}

# Export samples
step_output_summary <- function(ps, out_dir, type="unfiltered"){
  #Export raw csv
  phyloseq::psmelt(ps) %>%
    filter(Abundance > 0) %>%
    dplyr::select(-Sample) %>%
    write_csv(normalizePath(paste0(out_dir,"/raw_", type,".csv")))
  
  # Export species level summary of filtered results
  ps %>%
    phyloseq::psmelt() %>%
    filter(Abundance > 0) %>%
    left_join(refseq(ps) %>% as.character() %>% enframe(name="OTU", value="sequence")) %>%
    dplyr::select(OTU, sequence, rank_names(ps), sample_id, Abundance ) %>%
    pivot_wider(names_from = sample_id,
                values_from = Abundance,
                values_fill = list(Abundance = 0)) %>%
    write.csv(file = normalizePath(paste0(out_dir,"/summary_", type, ".csv")))

  #Output fasta of all ASV's
  seqs <- Biostrings::DNAStringSet(as.vector(phyloseq::refseq(ps)))
  Biostrings::writeXStringSet(seqs, filepath = normalizePath(paste0(out_dir,"/asvs_",type,".fasta")), width = 100) 

  out_paths <- c(
    paste0(out_dir,"/raw_", type,".csv"),
    paste0(out_dir,"/summary_",type,".csv"),
    paste0(out_dir,"/asvs_",type,".fasta")
  )
  if(!is.null(phy_tree(ps, errorIfNULL = FALSE))){
    #Output newick tree
    write.tree(phy_tree(ps), file= normalizePath(paste0(out_dir,"/tree_",type,".nwk")))
    out_paths <- c(out_paths, paste0(out_dir,"/tree_",type,".nwk"))
  }
  return(out_paths)
}

step_output_ps <- function(ps, out_dir, type="unfiltered"){
  seqtab <- phyloseq::otu_table(ps) %>%
    as("matrix") %>%
    as_tibble(rownames = "sample_id")
  
  taxtab <- phyloseq::tax_table(ps) %>%
    as("matrix") %>%
    as_tibble(rownames = "OTU") %>%
    unclassified_to_na(rownames = FALSE)
  
  #Check taxonomy table outputs
  if(!all(colnames(taxtab) == c("OTU", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))){
    message("Warning: Taxonomy table columns do not meet expectations for the staging database \n
          Database requires the columns: OTU, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species ")
  }
  
  samdf <- phyloseq::sample_data(ps) %>%
    as("matrix") %>%
    as_tibble()
  
  # Write out
  write_csv(seqtab, normalizePath(paste0(out_dir,"/seqtab_",type,".csv")))
  write_csv(taxtab, normalizePath(paste0(out_dir,"/taxtab_",type,".csv")))
  write_csv(samdf, normalizePath(paste0(out_dir,"/samdf_",type,".csv")))
  saveRDS(ps, paste0(out_dir,"/ps_",type,".rds"))
  
  out_paths <- c(
    paste0(out_dir,"/seqtab_",type,".csv"),
    paste0(out_dir,"/taxtab_",type,".csv"),
    paste0(out_dir,"/samdf_",type,".csv"),
    paste0(out_dir,"/ps_",type,".rds")
  )
  return(out_paths)
}



# Utilities ---------------------------------------------------------------
setup_multithread <- function(multithread, create_future = FALSE, quiet = FALSE) {
  ncores <- future::availableCores()
  if (isTRUE(multithread)) {
    cores <- ncores - 1
    if (!quiet) {
      message("Multithreading with ", cores, " cores")
    }
  } else if (is.numeric(multithread) & multithread > 1) {
    cores <- multithread
  } else if (isFALSE(multithread) | multithread == 1) {
    cores <- 1
  } else (stop("Multithread must be a logical or numeric vector of the numbers of cores to use"))
  
  if (cores > ncores) {
    cores <- ncores
    warning("The value provided to multithread is higher than the number of cores, using ", cores, " cores instead")
  } else {
    if (!quiet & cores > 1) {message("Multithreading with ", cores, " cores")}
  }
  
  # Handle multithread types
  if(create_future & cores >1){
    future::plan(future::multiprocess, workers = cores)
  } else if(create_future & cores == 1){
    future::plan(future::sequential)
  } 
  return(cores)
}

remove_if_exists <- function(file, quiet=FALSE){
  if(file.exists(file)) {
    if(file.remove(file)) {
      if(!quiet) message("Overwriting file:", file)
    } else {
      stop("Failed to overwrite file:", file)
    }
  }
}

create_samplesheet <- function(SampleSheet, runParameters, template = "V4"){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (length(SampleSheet) > 1) {multi <- TRUE}
  if (!length(SampleSheet) == length(runParameters)) {
    stop("Error: you have provided ", length(SampleSheet) , " SampleSheets and ", length(runParameters), " runParameters files. One of each must be provided per run")
  }
  
  #Parse files
  merged <- purrr::map2(SampleSheet, runParameters, parse_seqrun) %>%
    dplyr::bind_rows()
  
  # Reformat to the format required
  if (is.character(template) && template=="V4"){
    # Define template fields
    template_fields <- c("sample_id", "sample_name", "extraction_rep", "amp_rep", "client_name", 
                         "experiment_name", "sample_type", "collection_method", "collection_location", "lat_lon",
                         "environment", "collection_date", "operator_name", "description", "assay",
                         "extraction_method", "amp_method", "target_gene", "pcr_primers", "for_primer_seq",
                         "rev_primer_seq", "index_plate", "index_well", "i7_index_id", "i7_index",
                         "i5_index_id", "i5_index", "seq_platform", "fcid", "for_read_length", 
                         "rev_read_length", "seq_run_id", "seq_id", "seq_date", "analysis_method",
                         "notes"
    )
  } else if (any(class(template) == "data.frame")){
    template_fields <- colnames(template)
  } else {
    stop("Error, only template='V4' or a user provided data framecurrently supported")
  }
  
  # lookup table for renaming
  lookup <- c(i7_index = "index",
              i5_index = "index2",
              index_plate = "sample_plate",
              index_well = "sample_well",
              operator_name = "investigator_name",
              client_name = "project_name",
              seq_id = "instrument_name",
              seq_date = "run_start_date",
              seq_run_id = "run_id")
  
  matching <- merged %>%
    janitor::clean_names()%>%
    rename(any_of(lookup)) %>%
    dplyr::select_if(names(.) %in% template_fields)
  matching[,setdiff(template_fields, colnames(matching))] <- NA
  out <- matching %>%
    dplyr::select(all_of(template_fields))
  
  message(paste0(length(unique(out$sample_id))," samples total"))
  return(out)
}

parse_seqrun <- function(SampleSheet, runParameters){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (!length(SampleSheet) == length(runParameters)) {stop("Error: SampleSheet and RunParameters need to be provided for every run")}
  #detect format for run
  if(any(stringr::str_detect(readr::read_lines(runParameters), "MiSeq"))){
    format <- "miseq"
    sampleskip <- 20
    header_n_max <- 19
    reads_skip <- 12
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "novaseq"))){
    format <- "novaseq"
    sampleskip = 19
    header_n_max = 18
    reads_skip = 11
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "hiseq"))){
    format <- "hiseq"
    stop("Error: HiSeq not currently supported")
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "nextseq"))){
    format <- "nextseq"
    stop("Error: NextSeq not currently supported")
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "iseq"))){
    format <- "iseq"
    stop("Error: iSeq not currently supported")
  } else(
    stop("Error: compatable platfrom not detected in runParameters file")
  )
  # Read in samplesheet from run
  sample_sheet <- readr::read_csv(SampleSheet, skip=sampleskip, col_types = cols(
    Sample_ID = col_character(),
    Sample_Name = col_character(),
    Sample_Plate = col_character(),
    Sample_Well = col_character(),
    I7_Index_ID = col_character(),
    index = col_character(),
    I5_Index_ID = col_character(),
    index2 = col_character(),
    Sample_Project = col_character()
  ))
  
  withCallingHandlers({ # Handle Annoying missing columns function
    sample_header <- readr::read_csv(SampleSheet, n_max=header_n_max) %>%
      dplyr::select(1:2) %>%
      magrittr::set_colnames(c("var", "value")) %>%
      tidyr::drop_na(var) %>%
      dplyr::mutate(var = var %>%
                      str_replace_all("InvestigatorName", "Investigator_Name") %>% #Convert camel to snake case
                      str_replace_all("ExperimentName", "Experiment_Name") %>%
                      str_replace(" ", "_") %>%
                      make.unique() %>%
                      str_replace("\\.1", "_R")
      ) %>%
      tibble::column_to_rownames("var") %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::select_if(names(.) %in% c('Investigator_Name', 'Project_Name', 'Experiment_Name', 'Assay', 'Adapter'))
    
    reads <- readr::read_csv(SampleSheet, skip=reads_skip, n_max=2, col_types = cols_only(
      `[Reads]` = col_number() )) %>%
      pull(`[Reads]`)
    reads <- tibble::tibble(for_read_length = reads[1], rev_read_length = reads[2])
  },
  warning=function(w) {if (startsWith(conditionMessage(w), "Missing column names"))
    invokeRestart("muffleWarning")})
  
  # Read runparameters xml
  xmlFromRunParameters <- XML::xmlParse(runParameters)
  run_params <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "/RunParameters")) %>%
    as.data.frame(stringsAsFactors=FALSE )
  
  if(format == "miseq"){
    run_params <- run_params %>%
      dplyr::mutate(FlowCellExpiry = FlowcellRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    ReagentKitExpiry = ReagentKitRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    PR2Expiry = PR2BottleRFIDTag %>%
                      stringr::str_replace("^.{0,23}", "") %>%
                      stringr::str_replace(".{0,9}$", "") %>%
                      as.Date(),
                    FCID = Barcode %>%
                      stringr::str_replace("^.{0,10}", ""),
                    RunStartDate = lubridate::ymd(RunStartDate)
      ) %>%
      dplyr::rename( InstrumentName = ScannerID
      ) %>%
      dplyr::select(
        RunID,
        InstrumentName,
        RunNumber,
        FCID,
        RunStartDate,
        PR2BottleBarcode,
        ReagentKitBarcode,
        FlowCellExpiry,
        ReagentKitExpiry,
        PR2Expiry,
        MostRecentWashType) %>%
      dplyr::mutate_if(is.factor, as.character)
    
  } else if(format == "novaseq"){
    run_params <- run_params %>%
      dplyr::mutate(RunStartDate = lubridate::ymd(RunStartDate),
                    RunID = RunId
      ) %>%
      dplyr::select(
        RunID,
        InstrumentName,
        RunNumber,
        RunStartDate) %>%
      dplyr::mutate_if(is.factor, as.character)
    
    RFIDS <- XML::xmlToDataFrame(nodes = XML::getNodeSet(xmlFromRunParameters, "//RfidsInfo"))%>%
      as.data.frame(stringsAsFactors=FALSE) %>%
      dplyr::rename(
        FCID = FlowCellSerialBarcode,
        LibTubeID = LibraryTubeSerialBarcode,
        SbsID = SbsSerialBarcode,
        ClusterID = ClusterSerialBarcode,
        BufferID = BufferSerialBarcode,
      ) %>%
      dplyr::mutate(
        FlowCellExpiry = lubridate::mdy(str_remove(FlowCellExpirationdate," 00:00:00")),
        SbsExpiry = lubridate::mdy(str_remove(SbsExpirationdate," 00:00:00")),
        ClusterExpiry = lubridate::mdy(str_remove(ClusterExpirationdate," 00:00:00")),
        BufferExpiry = lubridate::mdy(str_remove(BufferExpirationdate," 00:00:00")),
      )%>%
      dplyr::select(
        FCID,
        LibTubeID,
        ClusterID,
        SbsID,
        BufferID,
        FlowCellExpiry,
        SbsExpiry,
        ClusterExpiry,
        BufferExpiry
      ) %>%
      dplyr::mutate_if(is.factor, as.character)
    run_params <- dplyr::bind_cols(run_params, RFIDS)
  }
  
  #Merge different sheets
  combined <- sample_sheet %>%
    cbind(sample_header) %>%
    cbind(reads) %>%
    cbind(run_params)
  message("Combined sample sheets for: ")
  message(paste0(unique(combined[]$FCID)," ", format, "\n"))
  return(combined)
}

download_zenodo <- function(doi, path = ".", quiet = FALSE) {

  # check for existence of the folder
  if(!dir.exists(path)){
    dir.create(path)
  }
  
  record <- str_remove(doi, fixed("10.5281/zenodo."))
  
  # Retrieve file name by records call
  base_url <- 'https://zenodo.org/api/records/'
  req <- curl::curl_fetch_memory(paste0(base_url, record))
  content <- jsonlite::fromJSON(rawToChar(req$content))
  
  # extract individual file names and urls
  file_urls <- content$files$links$self
  
  # extract check-sum(s)
  file_md5 <- content$files$checksum
  
  # Download files
  for (i in seq_along(file_urls)) {
    file_url <- file_urls[i]
    filename <- str_match(file_url, ".+/([^/]+)")[,2]
    destfile <- file.path(path, filename)
    
    # download file
    curl::curl_download(file_url, destfile, quiet=quiet)
    
    # Check file integrity
    md5 <- unname(tools::md5sum(destfile))
    zenodo_md5 <- str_split(file_md5[i], ":")[[1]][2]
    if (all.equal(md5, zenodo_md5)) {
      if (!quiet) message(filename," was downloaded and its integrity verified (md5sum: ", md5,")")
    } else {
      warning("Incorrect download! md5sum ", md5, " for file", filename, " does not match the Zenodo archived md5sum ", zenodo_md5)
    }
    
  }
}

new_bind <- function(a, b) {
  common_cols <- intersect(names(a), names(b))
  b[common_cols] <- map2_df(b[common_cols], 
                            map(a[common_cols], class), ~{class(.x) <- .y;.x})
  bind_rows(a, b)  
}

# Phyloseq utilities ------------------------------------------------------


phyloseq_filter_sample_wise_abund_trim <- function(physeq, minabund = 10, relabund = FALSE, rm_zero_OTUs = TRUE){
  
  ## Censore OTU abundance
  if(relabund == FALSE){     # trim based on absolute OTU counts
    
    res <- phyloseq::transform_sample_counts(physeq, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
    
  } else {                   # trim based on relative abundances within sample, but return original counts
    
    if(!minabund > 0 & minabund <= 1){
      stop("Error: for relative abundance trimmin 'minabund' should be in (0,1] interval.\n")
    }
    
    ## Convert data to relative abundances
    res <- phyloseq_standardize_otu_abundance(physeq, method = "total")
    
    ## Remove relative abundances less than the threshold value
    res <- phyloseq::transform_sample_counts(res, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
    
    ## Sample sums and data orientation
    smps <- phyloseq::sample_sums(physeq)
    if(phyloseq::taxa_are_rows(physeq) == TRUE){
      mar <- 2
    } else {
      mar <- 1
    }
    
    ## Convert back to counts by multiplying relative abundances by sample sums
    phyloseq::otu_table(res) <- phyloseq::otu_table(
      sweep(x = phyloseq::otu_table(res), MARGIN = mar, STATS = smps, FUN = `*`),
      taxa_are_rows = phyloseq::taxa_are_rows(physeq))
  }
  
  ## Remove zero-OTUs
  if(rm_zero_OTUs == TRUE){
    if (any(phyloseq::taxa_sums(res) > 0)){
      res <- phyloseq::prune_taxa(phyloseq::taxa_sums(res) > 0, res)
    } else {
      res <- NULL
    }
  }
  return(res)
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
    scale_x_continuous(labels =  label_number(scale_cut = cut_short_scale())) +
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