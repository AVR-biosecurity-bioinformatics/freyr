# Parse input arguments 
options <- commandArgs(trailingOnly = TRUE)
options[options ==""] <- NA_character_
options <- as.list(options)
#print(options)

names(options) <- c(
  "sample_sheet",
  "run_parameters",
  "read_dir",
  "pcr_primers",
  "for_primer_seq",
  "rev_primer_seq",
  "target_gene",
  "max_primer_mismatch",
  "read_min_length",
  "read_max_length",
  "read_max_ee",
  "read_trunc_length",
  "read_trim_left",
  "read_trim_right",
  "asv_min_length",
  "asv_max_length",
  "high_sensitivity",
  "concat_unmerged",
  "genetic_code",
  "coding",
  "phmm",
  "idtaxa_db",
  "ref_fasta",
  "idtaxa_confidence",
  "run_blast",
  "blast_min_identity",
  "blast_min_coverage",
  "target_kingdom",
  "target_phylum",
  "target_class",
  "target_order",
  "target_family",
  "target_genus",
  "target_species",
  "min_sample_reads",
  "min_taxa_reads",
  "min_taxa_ra"
)

print(options)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

# Define functions
create_samplesheet <- function(SampleSheet, runParameters, read_dirs, template = "V4"){
  if (missing(runParameters)) {stop("Error: need to provide a runParameters file in .xml format")}
  if (missing(SampleSheet)) {stop("Error: need to provide a SampleSheet file in .csv format")}
  if (length(SampleSheet) > 1) {multi <- TRUE}
  if (!length(SampleSheet) == length(runParameters)) {
    stop("Error: you have provided ", length(SampleSheet) , " SampleSheets and ", length(runParameters), " runParameters files. One of each must be provided per run")
  }
  
  # Parse all input samplesheets
  merged <- purrr::map2(SampleSheet, runParameters, parse_seqrun) %>%
    set_names(read_dirs) %>%
    bind_rows(.id="read_dir") %>%
    dplyr::bind_rows()
  
  # Reformat to internal samplesheet format
  if (is.character(template) && template=="V4"){
    # Define template fields
    template_fields <- c(
      "sample_id", "sample_name", "fcid", "target_gene", "pcr_primers", "for_primer_seq",
      "rev_primer_seq", "read_dir", 
      #"fwd", "rev", 
      "extraction_rep", "amp_rep", "client_name", 
      "experiment_name", "sample_type", "collection_method", "collection_location", "lattitude",
      "longitude", "environment", "collection_date", "operator_name", "description", "assay",
      "extraction_method", "amp_method", "index_plate", "index_well", "i7_index_id", "i7_index",
      "i5_index_id", "i5_index", "seq_platform", "for_read_length", "rev_read_length", 
      "seq_run_id", "seq_id", "seq_date", "analysis_method", "notes"
    )
  } else if (any(class(template) == "data.frame")){
    template_fields <- colnames(template)
  } else {
    stop("Error, only template='V4' or a user provided data framecurrently supported")
  }
  
  # lookup table for renaming columns
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
    dplyr::rename(any_of(lookup)) %>%
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
  
  #detect instrument used for run
  if(any(stringr::str_detect(readr::read_lines(runParameters), "MiSeq"))){
    format <- "miseq"
  } else if (any(stringr::str_detect(readr::read_lines(runParameters), "novaseq"))){
    format <- "novaseq"
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
  
  # Find positions of header and sample lines
  lines <- readr::read_lines(SampleSheet)
  #reads_line
  data_start_line <- which(str_detect(lines, "^\\[Data\\]"))
  header_end_line <- data_start_line-1
  reads_start_line <- which(str_detect(lines, "^\\[Reads\\]"))-1
  
  # Read in samplesheet from run
  sample_sheet <- readr::read_csv(SampleSheet, skip=data_start_line, col_types = cols(
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
  
  # Extract relevent info from sample header and add to sample sheet
  header_lines <- lines[1:header_end_line] %>%
    stringr::str_remove(",,.*$")
  
  sample_sheet$Investigator_Name <- if(any(str_detect(header_lines, "^Investigator Name"))){
    header_lines[which(str_detect(header_lines, "^Investigator Name"))] %>% str_remove("^.*,")
  } else {NA_character_}
  sample_sheet$Project_Name <- if(any(str_detect(header_lines, "^Project Name"))){
    header_lines[which(str_detect(header_lines, "^Project Name"))] %>% str_remove("^.*,")
  } else {NA_character_}  
  sample_sheet$Experiment_Name <- if(any(str_detect(header_lines, "^Experiment Name"))){
    header_lines[which(str_detect(header_lines, "^Experiment Name"))] %>% str_remove("^.*,")
  } else {NA_character_}  
  sample_sheet$Assay <- if(any(str_detect(header_lines, "^Assay"))){
    header_lines[which(str_detect(header_lines, "^Assay"))] %>% str_remove("^.*,")
  } else {NA_character_}  
  sample_sheet$Adapter <- if(any(str_detect(header_lines, "^Adapter"))){
    header_lines[which(str_detect(header_lines, "^Adapter"))] %>% str_remove("^.*,")
  } else {NA_character_} 
  sample_sheet$for_read_length <- if(any(str_detect(header_lines, "^\\[Reads\\]"))){
    header_lines[which(str_detect(header_lines, "^\\[Reads\\]"))+1] %>% str_remove("^.*,")
  } else {NA_character_}
  sample_sheet$rev_read_length <- if(any(str_detect(header_lines, "^\\[Reads\\]"))){
    header_lines[which(str_detect(header_lines, "^\\[Reads\\]"))+2] %>% str_remove("^.*,") 
  } else {NA_character_}  
  
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
  
  #Combine merge run params into sample sheet
  combined <- sample_sheet %>%
    cbind(run_params) 
  message("Combined sample sheets for: ")
  message(paste0(unique(combined[]$FCID)," ", format, "\n"))
  return(combined)
}

str_split_new <- function(string, pattern) {
  str_split(string, pattern)[[1]]
  }

# get number of primers
primer_length <- length(options$pcr_primers %>%
                          str_split_new(",|;"))

# Check that number of provided parameters are either 1, or match the number of provided primers
param_lengths <- options[!names(options) %in% c("sample_sheet", "run_parameters", "read_dir", "pcr_primers")] %>%
  purrr::map_dbl(~length(.x %>% str_split_new(",|;")))

if(any(!param_lengths %in% c(1, primer_length))){
  stop(names(param_lengths[param_lengths %in% c(1, primer_length)]), " parameters must have either 1 argument, or the same as pcr_primers (", primer_length,")")
}

# Get sample sheet & runparameters
sample_sheet <- options$sample_sheet%>%
  str_split_new(",") %>%
  normalizePath()
run_parameters <- options$run_parameters%>%
  str_split_new(",") %>%
  normalizePath()
read_dirs <- options$read_dir%>%
  str_split_new(",") %>%
  normalizePath()

message("Creating new params file from input arguments")

# Create parameters sheet
params <- options[!names(options) %in% c("sample_sheet", "run_parameters", "read_dir")] %>%
  purrr::map(~{
    #if(is.na(.x)){return(NA)}
    .x %>%
      str_split_new(",|;") %>%
      str_remove("\\[.*\\]")}) %>%
  bind_rows() %>%
  dplyr::select(any_of(c(
    "pcr_primers", "target_gene", "max_primer_mismatch", "read_min_length",
    "read_max_length", "read_max_ee", "read_trunc_length", "read_trim_left",
    "read_trim_right", "high_sensitivity", "asv_min_length", "asv_max_length",
    "concat_unmerged", "genetic_code", "coding", "phmm", "idtaxa_db",
    "ref_fasta", "idtaxa_confidence", "run_blast", "blast_min_identity",
    "blast_min_coverage", "target_kingdom", "target_phylum", "target_class",
    "target_order", "target_family", "target_genus", "target_species",
    "min_sample_reads", "min_taxa_reads", "min_taxa_ra")))

write_csv(params, "inputs/loci_params.csv")

# Create samplesheets
message("Creating new sample data file from input arguments")

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = sample_sheet, runParameters = run_parameters, read_dir=read_dirs, template = "V4") %>%
  distinct()

# Check that sample_ids contain fcid, if not; attatch
samdf <- samdf %>%
  mutate(sample_id = case_when(
    !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
    TRUE ~ sample_id
  ))

# Validate samples match the sample sheet
fastqFs <- purrr::map(read_dirs,
                      list.files, pattern="_R1_", full.names = TRUE) %>%
  unlist() %>%
  str_remove(pattern = "^(.*)\\/") %>%
  str_remove(pattern = "(?:.(?!_S))+$")

fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

missing_fastqFs <- setdiff(fastqFs, samdf$sample_id)
if (length(missing_fastqFs) > 0) warning("The fastq file/s: ", missing_fastqFs, " are not in the sample sheet")

missing_sample_ids <- setdiff(samdf$sample_id, fastqFs)
if (length(missing_sample_ids) > 0) {
  warning("The fastq file: ", missing_sample_ids, " is missing, dropping from samplesheet \n")
  samdf <- samdf %>% filter(!sample_id %in% missing_sample_ids)
}

# If there are brackets there, add primer details to the sample sheet based on string matching
pattern_matching <- any(str_detect(c(options$pcr_primers, options$for_primer_seq, options$rev_primer_seq), "\\[|\\]"))

if(pattern_matching){
  # Function to dynamically create case_when logic for primers
  create_case_when <- function(input) {
    case_when_expr <- map2_chr(names(input), input, ~ sprintf("str_detect(sample_id, '%s') ~ '%s'", .x, .y))
    case_when_expr <- paste(case_when_expr, collapse = ", ")
    case_when_expr <- paste0("case_when(", case_when_expr, ", TRUE ~ NA_character_)")
    rlang::parse_expr(case_when_expr)
  }
  
  # Helper function to create named vectors
  create_named_vector <- function(input) {
    input_values <- str_split_new(input, ",") %>% str_remove("\\[.*\\]")
    input_names <- str_split_new(input, ",") %>% str_remove("].*$") %>% str_remove("^.*\\[")
    setNames(input_values, input_names)
  }
  
  # Define named vectors for patterns and their corresponding primer sequences
  pattern_to_primers <- create_named_vector(options$pcr_primers)
  pattern_to_for <- create_named_vector(options$for_primer_seq)
  pattern_to_rev <- create_named_vector(options$rev_primer_seq)
  pattern_to_gene <- create_named_vector(options$target_gene)
  
  
  # Add primers to sample sheet based on sample_name matching
  samdf <- samdf %>%
    mutate(pcr_primers = !!create_case_when(pattern_to_primers),
           for_primer_seq = !!create_case_when(pattern_to_for),
           rev_primer_seq = !!create_case_when(pattern_to_rev),
           target_gene = !!create_case_when(pattern_to_gene))
  
} else if(!pattern_matching) {
  # Otherwise just add it direct 
  samdf <- samdf %>%
    mutate(
      pcr_primers = options$pcr_primers,
      for_primer_seq = options$for_primer_seq,
      rev_primer_seq = options$rev_primer_seq,
      target_gene = options$target_gene
    )
}


# Write out sample tracking sheet
write_csv(samdf, "inputs/Sample_info.csv")
