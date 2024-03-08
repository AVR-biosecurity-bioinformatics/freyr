#!/usr/bin/env Rscript

# run flow cell QC
step_seq_qc(flowcell_id)

step_switching_calc2 <- function(fcid, barcode_mismatch=1, quiet=FALSE){
  
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
  
  if (exists("gg.switch")) { 
    class(gg.switch)
    stop("gg.switch exists") 
    }

  pdf(file=paste(fcid,"_index_switching.pdf"), width = 11, height = 8 , paper="a4r")
    plot(gg.switch)
    try(dev.off(), silent=TRUE)
  return(res)
}




# run index switching calculation
step_switching_calc2(flowcell_id)

# copy output files to log folder
file.copy(paste0(flowcell_id,"_flowcell_qc.pdf"),paste0(projectDir,"/output/logs/",flowcell_id))
#file.copy(paste0(flowcell_id,"_index_switching.pdf"),paste0(projectDir,"/output/logs/",flowcell_id))


quit(status = 0)

# tar_target(seq_qc, {
#             process <- temp_samdf1_grouped %>%
#             dplyr::group_by(fcid) %>%
#             tidyr::nest() %>%
#             dplyr::mutate(seq_qc = purrr::map(fcid, step_seq_qc, quiet=FALSE, write_all=FALSE))
#             out <- paste0("output/logs/", unique(process$fcid),"/",unique(process$fcid),"_flowcell_qc.pdf")
#             if(is.na(process$seq_qc[[1]]$reads_pf)){
#             pdf(file=out, paper="A4")
#             plot.new()
#             text(x=.5, y=.5, "ERROR: InterOp folder or RunInfo.xml not present") 
#             try(dev.off(), silent=TRUE)
#             }
#             return(out)
#             }, pattern = map(temp_samdf1_grouped), format="file",  iteration = "vector"),

# tar_target(switching_qc,{
#             process <- temp_samdf1_grouped %>%
#             dplyr::group_by(fcid) %>%
#             tidyr::nest() %>%
#             dplyr::mutate(switching_qc = purrr::map(fcid, step_switching_calc, quiet=TRUE))
#             out <- paste0("output/logs/",unique(process$fcid),"/",unique(process$fcid),"_index_switching.pdf")
#             if(is.na(process$switching_qc[[1]]$switch_rate)){
#             pdf(file=out, paper="A4")
#             plot.new()
#             text(x=.5, y=.5, "ERROR: Undetermined reads file not found") 
#             try(dev.off(), silent=TRUE)
#             }
#             return(out)
#             },pattern = map(temp_samdf1_grouped), format="file", iteration = "vector"),
            