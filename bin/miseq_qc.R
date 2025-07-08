#!/usr/bin/env Rscript
tryCatch({

args <- R.utils::commandArgs(asValues = TRUE, trailingOnly = TRUE)

cat("\nArguments to process:\n")
str(args, no.list = T, nchar.max = 1E6)
cat("\n")

### process arguments 

read_group          <- args$read_group
miseq_dir           <- args$miseq_dir

sys.source(paste0(args$projectDir,"/bin/functions.R"), envir = .GlobalEnv)

### load only required packages
process_packages <- c(
    "dplyr",
    "ggplot2",
    "magrittr",
    "purrr",
    "readr",
    "savR",
    "seqateurs",
    "stringr",
    "tibble",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

#### TODO: Pull stats from top of index_switching.pdf and print to a final run report
### also flag index combos that have higher than average + (some stat)

#### TODO: use index_switch_calc.txt in jack_notes to run same process but in bash, which should be much faster

# extract "fcid" from read-header read group
fcid <- stringr::str_extract(read_group, "(?<=-)(\\S+)(?=__)", group = 1) 

# convert miseq_dir to a path
data_path <- base::normalizePath(miseq_dir)

seq_dir <- paste0(data_path,"/",fcid)
  
barcode_mismatch <- 1

# Check that required files exist
if(!dir.exists(seq_dir)) {
    stop("seq_dir doesnt exist, check that the correct path was provided")
}

if(!dir.exists(paste0(seq_dir,"/InterOp"))){
    stop("InterOp folder must be present to run quality checks")
}

if(!file.exists(paste0(seq_dir,"/RunInfo.xml"))){
    stop("RunInfo.xml must be present to run quality checks")
}

### run flow cell QC

# Sequencing run quality check using savR
fc <- savR::savR(seq_dir)

# Ensure indices are present
if(length(fc@parsedData)==0){
    stop(paste0("Flow cell metrics could not be parsed for", read_group))
}

# write output tables
readr::write_csv(savR::correctedIntensities(fc), paste0("correctedIntensities_",read_group,".csv"))
readr::write_csv(savR::errorMetrics(fc), paste0("errorMetrics_",read_group,".csv"))
readr::write_csv(savR::extractionMetrics(fc), paste0("extractionMetrics_",read_group,".csv"))
readr::write_csv(savR::qualityMetrics(fc), paste0("qualityMetrics_",read_group,".csv"))
readr::write_csv(savR::tileMetrics(fc), paste0("tileMetrics_",read_group,".csv"))

# plot flowcell qc
gg.avg_intensity <- 
    fc@parsedData[["savCorrectedIntensityFormat"]]@data %>%
    dplyr::group_by(tile, lane) %>%
    dplyr::summarise(Average_intensity = mean(avg_intensity), .groups="drop") %>% 
    dplyr::mutate(side = case_when(
        stringr::str_detect(tile, "^11") ~ "Top",
        stringr::str_detect(tile, "^21") ~ "Bottom"
    )) %>%
    ggplot2::ggplot(aes(x=lane, y=as.factor(tile), fill=Average_intensity)) +
    geom_tile() +
    facet_wrap(~side, scales="free") +
    scale_fill_viridis_c()

pdf(file=paste0("flowcell_qc_",read_group,".pdf"), width = 11, height = 8 , paper="a4r")
plot(gg.avg_intensity)
savR::pfBoxplot(fc)
for (lane in 1:fc@layout@lanecount) {
    savR::qualityHeatmap(fc, lane, 1:fc@directions)
}
try(dev.off(), silent=TRUE)

# Return the total number of reads and passing filter
readspassing <- 
    fc@parsedData[["savTileFormat"]]@data %>%
    dplyr::filter(code %in% c(100,101)) %>%
    dplyr::mutate(code = case_when(
        code == 100 ~ "reads_total",
        code == 101 ~ "reads_pf"
    ),read_group = stringr::str_remove(fc@flowcell, "^.*-")) %>% 
    dplyr::group_by(read_group, code) %>%
    dplyr::summarise(reads = sum(value), .groups="drop") %>%
    tidyr::pivot_wider(names_from = code,
                        values_from = reads)

readr::write_csv(readspassing, paste0("readsPassing_",read_group,".csv"))



### run index switching calculation

# Check if undetermined reads file exists
if(!any(stringr::str_detect(list.files(seq_dir, pattern="_R1_", full.names = TRUE), "Undetermined"), na.rm = TRUE)){
    stop("Error, an Undetermined reads fastq must be present to calculate index switching")
}

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
    stop(paste0("No index sequences present in fastq headers for run", read_group, " no switch rate calculated"))
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
    dplyr::group_modify(~{
    .x %>%
        dplyr::top_n(n=1, Freq) %>%
        dplyr::slice(1) %>%  # Handle ties
        dplyr::mutate(Freq = sum(.x$Freq)) 
    })

# Check if indices are combinatorial
if(any(duplicated(applied_indices$index)) | any(duplicated(applied_indices$index2))){
    stop(paste0("Combinatorial indexes detected for", read_group, " no switch rate calculated"))
}

# Get other undetermined reads which had completely unapplied indexes
other_reads <- dplyr::anti_join(indices,combos, by=c("index", "index2")) %>%
    dplyr::summarise(sum = sum(Freq, na.rm = TRUE)) %>%
    dplyr::pull(sum)

#Summary of index switching rate
res <- 
    switched %>%
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

message("Index switching rate calculated as: ", res$switch_rate)

#Plot switching - handling barcode mismatch during demultiplexing 

# Update indexes using their hamming distance to those originally appliex
switch_plot_dat <- 
    switched %>%
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
    index_dist <- stringdist::stringdist(.x,index_list) # stringdist is imported by seqateurs
    # Remove those above barcode_mismatch threshold
    index_list <- index_list[index_dist <= barcode_mismatch]
    index_dist <- index_dist[index_dist <= barcode_mismatch]
    return(index_list[which.min(index_dist)])
    })) %>%
    tidyr::unnest(c(index, index2)) %>%
    dplyr::group_by(index, index2, Sample_Name) %>%
    dplyr::summarise(Freq = sum(Freq))

gg.switch <- 
    switch_plot_dat %>%
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
    labs(title= read_group, subtitle = paste0(
    "Total Reads: ", sum(indices$Freq, na.rm=TRUE),
    ", Switch rate: ", sprintf("%1.4f%%", res$switch_rate*100),
    ", Contam rate: ", sprintf("%1.6f%%", res$contam_rate*100),
    ", Other reads: ", other_reads)) 

pdf(file=paste0("index_switching_",read_group,".pdf"), width = 11, height = 8 , paper="a4r")
    plot(gg.switch)
    try(dev.off(), silent=TRUE)


# stop(" *** stopped manually *** ") ##########################################
}, 
finally = {
    ### save R environment if script throws error code
    if (args$rdata == "true") {save.image(file = paste0(args$process_name,".rda"))}
})