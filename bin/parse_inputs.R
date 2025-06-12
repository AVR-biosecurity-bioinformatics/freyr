#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "dplyr",
    "magrittr",
    "purrr",
    "readr",
    "stringr",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "launchDir",
    "projectDir",
    "params_dict",
    "samplesheet",
    "primer_params",
    "seq_type",
    "paired",
    "subsample"
)
lapply(nf_vars, nf_var_check)

### process variables 

### run code

# read in samplesheet tibble
samplesheet_df <- 
    readr::read_csv(samplesheet, show_col_types = F, skip_empty_rows = TRUE, na = c("", "NA", " ")) %>%
    dplyr::select(!tidyselect::starts_with("...")) %>% # remove any blank columns that have been named "...X"
    dplyr::filter(dplyr::if_any(tidyselect::everything(), ~ !is.na(.))) # remove completely blank lines

# read in primer parameters tibble
primer_params_df <- 
    readr::read_csv(primer_params, show_col_types = F, skip_empty_rows = TRUE, na = c("", "NA", " ")) %>%
    dplyr::select(!tidyselect::starts_with("...")) %>% # remove any blank columns that have been named "...X"
    dplyr::filter(dplyr::if_any(tidyselect::everything(), ~ !is.na(.))) # remove completely blank lines

### validation of samplesheet content

sse <- "SAMPLESHEET FILE ERROR: "
prpe <- "PRIMER PARAMS FILE ERROR: "
ppe <- "PIPELINE PARAMETER ERROR: "

### validate samplesheet contents ----------------------------------------------------------------------------------------


ssc <- colnames(samplesheet_df)

# check required columns are present
if( !any(c("sample", "primers") %in% ssc) ){
    stop(paste0(sse, "Samplesheet must contain both of the fields 'sample' and 'primers'!"))
}

read_col_which <- c("read_dir","fwd", "rev", "single")[which(c("read_dir","fwd", "rev", "single") %in% ssc)]

if( ! (read_col_which == "read_dir" || setequal(read_col_which ,c("fwd","rev")) || read_col_which == "single" ) ){
    stop(paste0(sse, "Samplesheet must contain only one of the following fields/field groups: 'read_dir'; 'fwd' AND 'rev'; or 'single'."))
}

# check if disallowed 'paired' and read column values
if ( paired == "true" && read_col_which == "single" ){
    stop(paste0(sse, "'single' column cannot be used with parameter '--paired'."))
}

if ( paired == "false" && setequal(read_col_which, c("fwd","rev"))){
    stop(paste0(sse, "'fwd' and 'rev' cannot be used with '--paired false'."))
}

# check sample field is only unique values
if (any(duplicated(samplesheet_df$sample))) { 
    sample_duplicates <- samplesheet_df$sample[which(duplicated(samplesheet_df$sample))]
    stop(paste0(sse, "The following 'sample' values are repeated: ",stringr::str_flatten(sample_duplicates, collapse = ", ") ,"\nAll 'sample' values in samplesheet must be unique."))
} 

## any spaces or commas in fields?
# if 'read_group' is present, validate it too, otherwise just validate other columns
if ("read_group" %in% ssc){
    if ( any(c(samplesheet_df$sample, samplesheet_df$primers, samplesheet$read_group) %>% stringr::str_detect(., "^[^\\s,]+$", negate = T)) ){
        stop(paste0(sse, "'sample', 'primers' and 'read_group' values must not contain spaces or commas."))
    } 
} else {
    if ( any(c(samplesheet_df$sample, samplesheet_df$primers) %>% stringr::str_detect(., "^[^\\s,]+$", negate = T)) ){
        stop(paste0(sse, "'sample' and 'primers'values must not contain spaces or commas."))
    } 
}

# validate read path fields (without validating reads themselves) -- no spaces or commas
if ( read_col_which == "read_dir" ){
    if (any(samplesheet_df$read_dir %>% stringr::str_detect(., "^[^\\s,]+$", negate = T))){
        stop(paste0(sse, "'read_dir' values must not contains spaces or commas."))
    }
} else if ( setequal(read_col_which ,c("fwd","rev")) ){
    if (any(samplesheet_df$fwd %>% stringr::str_detect(., "^[^\\s,]+$", negate = T))){
        stop(paste0(sse, "'fwd' values must not contains spaces or commas."))
    }
    if (any(samplesheet_df$rev %>% stringr::str_detect(., "^[^\\s,]+$", negate = T))){
        stop(paste0(sse, "'rev' values must not contains spaces or commas."))
    }
} else if ( read_col_which == "single" ){
    if (any(samplesheet_df$single %>% stringr::str_detect(., "^[^\\s,]+$", negate = T))){
        stop(paste0(sse, "'single' values must not contains spaces or commas."))
    }
} else {
    stop(paste0(sse, "Invalid read field"))
}


### parse and/or detect read paths
if ( read_col_which == "read_dir" & paired == "true" ) { # if "read_dir" column is used and reads are paired
    ## try to find reads
    # convert read directory to absolute path
    samplesheet_df_rp <- 
        samplesheet_df %>% 
        dplyr::mutate(
            read_dir = dplyr::case_when(
                stringr::str_starts(read_dir, "/") ~ read_dir, # if already absolute path, leave it be
                stringr::str_starts(read_dir, "\\./") ~ stringr::str_replace(read_dir, "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(read_dir, "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    reads_list = list() # create empty list
    for (i in 1:length(samplesheet_df_rp$read_dir)) { # loop through rows of samplesheet
        i_readfiles <- list.files( # find full paths of files matching sample with a fastQ extension
            path = samplesheet_df_rp$read_dir[i],
            pattern = paste0(samplesheet_df_rp$sample[i],"(?:_[^\\s/]+)?\\.f(ast)?q(\\.gz)?$"),
            full.names = T, 
            recursive = T
            ) %>% unlist()
        # check exactly two read files are found
        if ( length(i_readfiles) != 2 ) { 
            stop (paste0(sse, "Found ",length(i_readfiles)," read files matching '",samplesheet_df_rp$sample[i],"' in '",samplesheet_df_rp$read_dir[i],"' when 2 were expected.\n\tCheck you have filled out samplesheet correctly."))  
        }
        if (params.extension != "null") { # using params.extension if supplied
            message(paste0("Using 'params.extension' (",params.extension,") to find read files in supplied directories ('read_dir')."))
            # check read files match params.extension
            if (any(stringr::str_detect(i_readfiles, pattern = paste0(params.extension,"$")))) { 
                stop(paste0(sse, "Read files found in 'read_dir' do not match 'params.extension' file extension--check samplesheet and pipeline parameters."))
            }
            # sort read files by text after extension
            ### TODO: do this
        } else { # not using params.extension
            # sort read files by natural order (ie. 1 before 2; F before R)
            i_readfiles <- stringr::str_sort(i_readfiles)
        }
        
        reads_list[[i]] <- c(samplesheet_df_rp$sample[i], i_readfiles) # append to list
    }
    reads_df <- do.call(rbind, reads_list) # combine list into dataframe
    colnames(reads_df) <- c("sample","fwd","rev") # name columns
    reads_df <- 
        reads_df %>% tibble::as_tibble() # convert to tibble
    # add read paths to samplesheet
    samplesheet_df_rp <- 
        samplesheet_df_rp %>% dplyr::left_join(., reads_df, by = "sample") 

} else if ( setequal(read_col_which, c("fwd","rev"))) { # if 'fwd' and 'rev' columns both exist...
    # convert read paths to absolute paths
    samplesheet_df_rp <- samplesheet_df %>% 
        dplyr::mutate(
            across(
                c(fwd, rev),
                ~ dplyr::case_when(
                    stringr::str_starts(., "/") ~ ., # if already absolute path, leave it be
                    stringr::str_starts(., "\\./") ~ stringr::str_replace(., "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                    .default = stringr::str_replace(., "^", paste0(launchDir,"/")) # else append launchDir to front to produce absolute path
                )
            )
        )
} else if ( read_col_which == "read_dir" & paired == "false" ) { # if "read_dir" column exists in samplesheet
    ## try to find reads
    # convert read directory to absolute path
    samplesheet_df_rp <- samplesheet_df %>% 
        dplyr::mutate(
            read_dir = dplyr::case_when(
                stringr::str_starts(read_dir, "/") ~ read_dir, # if already absolute path, leave it be
                stringr::str_starts(read_dir, "\\./") ~ stringr::str_replace(read_dir, "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(read_dir, "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    reads_list = list() # create empty list
    for (i in 1:length(samplesheet_df_rp$read_dir)) { # loop through rows of samplesheet
        i_readfiles <- list.files( # find full paths of files matching sample with a fastQ extension
            path = samplesheet_df_rp$read_dir[i],
            pattern = paste0(samplesheet_df_rp$sample[i],"(?:_[^\\s/]+)?\\.f(ast)?q(\\.gz)?$"),
            full.names = T, 
            recursive = T
        ) %>% unlist()
        # check exactly one read file is found per sample
        if ( length(i_readfiles) != 1 ) { 
            stop (paste0(sse, "Found ",length(i_readfiles)," read files matching '",samplesheet_df_rp$sample[i],"' in '",samplesheet_df_rp$read_dir[i],"' when 1 was expected.\n\tCheck you have filled out samplesheet correctly."))  
        }
        if (params.extension != "null") { # using params.extension if supplied
            message(paste0("Using 'params.extension' (",params.extension,") to find read files in supplied directories ('read_dir')."))
            # check read files match params.extension
            if (any(stringr::str_detect(i_readfiles, pattern = paste0(params.extension,"$")))) { 
                stop(paste0(sse, "Read files found in 'read_dir' do not match 'params.extension' file extension--check samplesheet and pipeline parameters."))
            }
            # sort read files by text after extension
            ### TODO: do this
        } else { # not using params.extension
            # sort read files by natural order (ie. 1 before 2; F before R)
            i_readfiles <- stringr::str_sort(i_readfiles)
        }
        
        reads_list[[i]] <- c(samplesheet_df_rp$sample[i], i_readfiles) # append to list
    }
    reads_df <- do.call(rbind, reads_list) # combine list into dataframe
    colnames(reads_df) <- c("sample","single") # name columns
    reads_df <- 
        reads_df %>% tibble::as_tibble() # convert to tibble
    # add read paths to samplesheet
    samplesheet_df_rp <- 
        samplesheet_df_rp %>% dplyr::left_join(., reads_df, by = "sample") 

} else if (( read_col_which == "single" )) { # if single reads only
# convert read paths to absolute paths
    samplesheet_df_rp <- samplesheet_df %>% 
        dplyr::mutate(
            dplyr::across(
                c(single),
                ~ dplyr::case_when(
                    stringr::str_starts(., "/") ~ ., # if already absolute path, leave it be
                    stringr::str_starts(., "\\./") ~ stringr::str_replace(., "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                    .default = stringr::str_replace(., "^", paste0(launchDir,"/")) # else append launchDir to front to produce absolute path
                )
            )
        )
} else {
    stop (paste0(sse, "Invalid samplesheet read columns and 'paired' pipeline parameter."))
}

# check that read file paths are unique
if (paired == "true") {
    if (any(duplicated(samplesheet_df_rp$fwd))) {stop (paste0(sse, "At least two samples share the same forward read file in the samplesheet!"))}
    if (any(duplicated(samplesheet_df_rp$rev))) {stop (paste0(sse, "At least two samples share the same reverse read file in the samplesheet!"))}
} else if ( paired == "false" ) {
    if (any(duplicated(samplesheet_df_rp$single))) {stop (paste0(sse, "At least two samples share the same read file in the samplesheet!"))}
} else {
    stop(paste0(ppe, "'paired' = '", paired, "' but must be 'true' or 'false'."))
}

### extract read_group from read headers

# # example for MiSeq data
# id_miseq <- "M01054:780:000000000-K77JP:1:2116:15731:20801 1:N:0:TGGCATGT+CCTTGTAG"
# # novaseq
# id_nova <- "A00878:92:HTFL5DRX2:2:2101:9272:1016 2:N:0:ACACGGGCCA+GTTGACCGTT"
# # internal (converted) MGI example 
# id_mgi_int <- "M:0:V350254379:4:C001R012:0:2065 1:N:0"
# # raw MGI example
# id_mgi_raw <- "V350231845L1C001R00100001249/1"
# # above interpretable as: "M:0:V350231845:1:1249:1: 1:N:0"
# # MinION format
# id_minion <- "837ed09c-2662-426a-bac9-ba59af484392 runid=3639b6d214857c63e9ea16dcba8f6471fbee372c ch=365 start_time=2025-03-05T14:38:32.214861+11:00 flow_cell_id=FBA86592 basecall_gpu=NVIDIA_RTX_A2000_12GB protocol_group_id=20250305_Grains_COI sample_id= barcode=barcode01 barcode_alias=barcode01 parent_read_id=837ed09c-2662-426a-bac9-ba59af484392 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v4.3.0"


### extract read_group from read headers if not present
if ("read_group" %in% colnames(samplesheet_df_rp)){
  
  if (any(samplesheet_df_rp$read_group %>% is.na())){
        stop(paste0("Some values of 'read_group' are NA!"))
  } else {
        # keep as-is -- information given by user
        samplesheet_df_rg <- samplesheet_df_rp
  }

} else {
  
    # extract read_group from read headers
    samplesheet_df_rg <- samplesheet_df_rp %>% dplyr::mutate(read_group = NA)

    ## loop through each sample
    for ( i in 1:length(samplesheet_df_rg$sample)){
        sample_i_name <- samplesheet_df_rg$sample[i]
        if ( paired == "true"){
            # if paired, use fwd read
            read_path <- samplesheet_df_rg$fwd[i]
        } else if ( paired == "false" ){
            # if unpaired, use single read
            read_path <- samplesheet_df_rg$single[i]
        }
        # sample 10 reads from read file
        read_sampler <- ShortRead::FastqSampler(read_path, 10)
        set.seed(1); read_yield <- ShortRead::yield(read_sampler)
        close(read_sampler) # close connection to file
        # get vector of headers
        read_headers <- read_yield@id %>% as.character
        
        ## detect format of sequencer and parse header accordingly
        # Illumina format (MiSeq or NovaSeq)
        if ( all(stringr::str_detect(read_headers, "^[[:alnum:]]+:\\d+:[[:alnum:]-]+:\\d+:\\d+:\\d+:\\d+ [12]:[YN]:\\d+(:[[:alpha:]+]+)?$"))){
        
            message(paste0("Read header format for sample '",sample_i_name,"' detected as Illumina (MiSeq or NovaSeq)"))
            if ( seq_type == "nanopore" ){
                stop(paste0(ppe, "'--seq_type' given as 'nanopore', but read header format appears to be Illumina!"))
            }
            # get sample flowcell (3rd element from : delimited list)
            sample_flowcell <- read_headers %>% stringr::str_split(., ":") %>% sapply(., `[`, 3) %>% unique()
            # check exactly one flowcell ID found
            if (length(sample_flowcell) != 1){
                stop(paste0("Invalid flowcell detection for sample '",sample_i_name,"'"))
            }
            # get sample lane (4th element from : delimited list)
            sample_lane <- read_headers %>% stringr::str_split(., ":") %>% sapply(., `[`, 4) %>% unique()
            # check exactly one lane number found
            if (length(sample_lane) != 1){
                stop(paste0("Invalid lane detection for sample '",sample_i_name,"'"))
            }
            # combine flowcell and lane into a single 'read_group' string
            sample_rg <- paste0(sample_flowcell,":",sample_lane)
        
        # MGI internal format 
        } else if ( all(stringr::str_detect(read_headers, "^[[:alnum:]]+:\\d+:[[:alnum:]-]+:\\d+:[[:alnum:]]+:\\d+:\\d+ [12]:[YN]:\\d+(:[[:alpha:]+]+)?$")) ) {
        
            message(paste0("Read header format for sample '",sample_i_name,"' detected as MGI (converted)"))
            if ( seq_type == "nanopore" ){
                stop(paste0(ppe, "'--seq_type' given as 'nanopore', but read header format appears to be MGI!"))
            }
            # get sample flowcell (3rd element from : delimited list)
            sample_flowcell <- read_headers %>% stringr::str_split(., ":") %>% sapply(., `[`, 3) %>% unique()
            # check exactly one flowcell ID found
            if (length(sample_flowcell) != 1){
                stop(paste0("Invalid flowcell detection for sample '",sample_i_name,"'"))
            }
            # get sample lane (4th element from : delimited list)
            sample_lane <- read_headers %>% stringr::str_split(., ":") %>% sapply(., `[`, 4) %>% unique()
            # check exactly one lane number found
            if (length(sample_lane) != 1){
                stop(paste0("Invalid lane detection for sample '",sample_i_name,"'"))
            }
            # combine flowcell and lane into a single 'read_group' string
            sample_rg <- paste0(sample_flowcell,":",sample_lane)
            
        # MGI raw format
        } else if ( all(stringr::str_detect(read_headers, "^[[:alnum:]]+/[12]$")) ) {
        
            message(paste0("Read header format for sample '",sample_i_name,"' detected as MGI (native)"))
            if ( seq_type == "nanopore" ){
                stop(paste0(ppe, "'--seq_type' given as 'nanopore', but read header format appears to be MGI!"))
            }
            # get sample flowcell (3rd element from : delimited list)
            sample_flowcell <- stringr::str_extract(read_headers, "^[[:alnum:]]+?(?=L)") %>% unique()
            # check exactly one flowcell ID found
            if (length(sample_flowcell) != 1){
                stop(paste0("Invalid flowcell detection for sample '",sample_i_name,"'"))
            }
            # get sample lane (4th element from : delimited list)
            sample_lane <- stringr::str_extract(read_headers, "(?<=L)\\d+") %>% unique()
            # check exactly one lane number found
            if (length(sample_lane) != 1){
                stop(paste0("Invalid lane detection for sample '",sample_i_name,"'"))
            }
            # combine flowcell and lane into a single 'read_group' string
            sample_rg <- paste0(sample_flowcell,":",sample_lane)
        
        # ONT format
        } else if ( all(stringr::str_detect(read_headers, " runid=[[:alnum:]]+")) && all(stringr::str_detect(read_headers, " flow_cell_id=[[:alnum:]]+"))  ) {
        
            message(paste0("Read header format for sample '",sample_i_name,"' detected as ONT"))
            if ( seq_type != "nanopore" ){
                stop(paste0(ppe, "'--seq_type' given as '",seq_type,"', but read header format appears to be ONT/'nanopore'!"))
            }

            # extract flowcell ID as read group (no lanes)
            sample_rg <- stringr::str_extract(id_minion, "(?<=flow_cell_id=)[[:alnum:]]+(?=( |$))") %>% unique()
            # check only one sample_rg extracted
            if (length(sample_rg) != 1){
                stop(paste0("Invalid flowcell detection for sample '",sample_i_name,"'"))
            }
        
        } else {
            stop(paste0("Read header format for sample '",sample_i_name,"' could not be parsed!"))
        }
        
        # update samplesheet tibble with new read_group for this sample
        samplesheet_df_rg <- 
            samplesheet_df_rg %>%
            dplyr::mutate(read_group = dplyr::if_else(sample == sample_i_name, sample_rg, read_group))
        
    }

}

# check if any 'read_group' values are still NA
if (any(samplesheet_df_rg$read_group %>% is.na())){
    stop(paste0(sse, "Not all 'read_group' values have been determined!"))
}

# pseudorandomly subsample samplesheet if params.subsample is defined
# this is done per pcr_primer x read_group combination so expected combinations are (likely) retained
if (subsample != "false") {
    { set.seed(1); samplesheet_df_ss <- 
        samplesheet_df_rg %>% 
        dplyr::group_by(primers, read_group) %>% 
        dplyr::slice_sample(n = as.numeric(subsample), replace = F) %>%
        dplyr::ungroup() }
} else {
    samplesheet_df_ss <- samplesheet_df_rg
}

# enforce column order
if ( paired == "true" ){
    samplesheet_complete <- 
        samplesheet_df_ss %>%
        dplyr::relocate(sample, read_group, primers, read_dir, fwd, rev)
} else {
    samplesheet_complete <- 
        samplesheet_df_ss %>%
        dplyr::relocate(sample, read_group, primers, read_dir, single)
}

## create samplesheet outputs
# create arbitrary metadata samplesheet -- keep only sample and read_group
samplesheet_complete %>%
    dplyr::select(-c(primers, read_dir, tidyselect::any_of(c("read_dir", "fwd", "rev", "single")))) %>%
    readr::write_csv(., "sample_metadata.csv")

# create samplesheet without 'sample_primers' field, for FASTQC channel
samplesheet_complete %>%
    dplyr::select(c(sample, read_group, primers, tidyselect::any_of(c("fwd","rev", "single")))) %>%
    readr::write_csv(., "samplesheet_unsplit.csv")

# split samplesheet by primers to generate samplesheet for split samples (without arbitrary metadata) -- main sample channel data
samplesheet_complete %>%
    dplyr::select(c(sample, read_group, primers, tidyselect::any_of(c("fwd","rev", "single")))) %>%
    tidyr::separate_longer_delim(primers, delim = ";") %>%
    dplyr::mutate(sample_primers = paste0(sample,"_",primers), .after = primers) %>%
    readr::write_csv(., "samplesheet_split.csv")



### validate primer_params content -------------------------------------------------------------------------------------

# check primers is correctly formatted
if (!all(stringr::str_detect(primer_params_df$primers, "^((?!__)[^\\s,;])+$"))){
    stop(paste0(prpe, "'primers' values must not contain spaces, commas, semi-colons or two underscores in a row."))
}

# check none of the primer_params values contain spaces or commas 
if (!all(stringr::str_detect(primer_params_df %>% unlist %>% unname, "^[^\\s,;]+$") %>% .[!is.na(.)] )){
    stop(paste0(prpe, "Primer parameters must not contain spaces, commas or semi-colons."))
}

# check primers column contains only unique values
if (any(duplicated(primer_params_df$primers))) {
    stop (paste0(prpe, "'primers' values are not unique!"))
}

# add missing columns with default values
default_params <- 
    tibble::tibble(
        primers = NA_character_,
        locus = NA_character_,
        for_primer_seq = NA_character_,
        rev_primer_seq = NA_character_,
        max_primer_mismatch = 0,
        read_min_length = 20,
        read_max_length = Inf,
        read_max_ee = 1,
        read_trunc_length = 0,
        read_trim_left = 0,
        read_trim_right = 0,
        asv_min_length = 0,
        asv_max_length = Inf,
        concat_unmerged = FALSE,
        genetic_code = NA_character_,
        coding = FALSE,
        phmm = NA_character_,
        idtaxa_db = NA_character_,
        ref_fasta = NA_character_,
        idtaxa_confidence = 60,
        run_blast = FALSE,
        blast_min_identity = 97,
        blast_min_coverage = 90,
        target_kingdom = NA_character_,
        target_phylum = NA_character_,
        target_class = NA_character_,
        target_order = NA_character_,
        target_family = NA_character_,
        target_genus = NA_character_,
        target_species = NA_character_,
        min_sample_reads = 0,
        min_taxa_reads = 0,
        min_taxa_ra = 0
    )

# combine default values and supplied values
primer_params_da <- 
    primer_params_df %>% 
    # add missing columns with default values (https://stackoverflow.com/a/53691922)
    tibble::add_column(
        ., 
        !!!default_params[setdiff(names(default_params), names(.))]
    ) %>%
    # reorder columns
    dplyr::relocate(tidyselect::all_of(names(default_params))) %>%
    # make every column characters, will convert later
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.character))


## check primers values match between samplesheet and primer_params
# get vector of unique primers values, undelimited from samplesheet
ss_pp <- samplesheet_complete$primers %>% stringr::str_split(., ";") %>% unlist() %>% unique()

if ( !setequal(primer_params_da$primers, ss_pp)){
    stop(paste0("SAMPLESHEET/PRIMER PARAMS ERROR: 'primers' values in samplesheet ('",stringr::str_flatten(ss_pp, " / "),"') and primer_params ('",stringr::str_flatten(primer_params_da$primers, " / "),"') do not match."))
}

### replace primer params from .csv with those from the '--pp_*' pipeline parameters
# vector of primer parameters
pp_vec <- 
    c(
        "max_primer_mismatch",
        "read_min_length",
        "read_max_length",
        "read_max_ee",
        "read_trunc_length",
        "read_trim_left",
        "read_trim_right",
        "asv_min_length",
        "asv_max_length",
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
    ) %>%
    stringr::str_replace(., "^", "pp_")

# use params_list from process_start.R to create a named vector of parameters and values
params_vec <- 
    sapply(
        params_list,
        function(x){
            y <- x[2]
            names(y) <- x[1]
            return(y)
        }
    )

# just the pp_ parameters 
params_pp <- params_vec[names(params_vec) %in% pp_vec]

# remove parameters set to "null" (ie. unset)
params_pp_nn <- 
    params_pp %>% stringr::str_subset(., "^null$", negate = TRUE)

# get vector of primers values
primers_vec <- primer_params_da$primers

# create new version of primer_params that will be updated through the loop
primer_params_up <- primer_params_da

if ( length(params_pp_nn) > 0 ){
    ## loop through non-null parameters updating the primer_params tibble
    for ( i in 1:length(params_pp_nn)) {
        
        # get name of parameter
        p_name <- names(params_pp_nn[i]) %>% stringr::str_remove(., "^pp_")
        # get value of parameter
        p_value <- params_pp_nn[i] %>% unname()
        
        # check format and make changes to data frame accordingly 
        if ( stringr::str_detect(p_value, "^[^\\[\\]\\;]+$")  ) { # 1 -- one value with no primers tag 
        
            # replace parameter with new value for all rows
            primer_params_up <- 
                primer_params_up %>% dplyr::mutate("{p_name}" := p_value)
            
            message(paste0("Primer parameter '",p_name,"' has been updated to '",p_value,"'."))
        
        } else if ( stringr::str_detect(p_value, "^\\[[^\\[\\]\\;]+\\][^\\[\\]\\;]+$") ) { # [A]1 -- one value with a primers tag 
        
            # get tag and non-tag parts of the parameter value
            p_tag <- stringr::str_extract(p_value, "^\\[(\\S+)\\]", group = 1)
            p_notag <- stringr::str_extract(p_value, "^\\[\\S+\\](\\S+)", group = 1)
            
            # check tag is a valid primers value
            if ( !p_tag %in% primers_vec ){
                stop(paste0(ppe, "The value of parameter '--pp_",p_name,"' does not have a valid primers tag: '", p_value,"'"))
            }
            
            # replace parameter just for given primers value
            primer_params_up <- 
                primer_params_up %>%
                dplyr::mutate(
                    dplyr::across(
                        {{p_name}},
                        ~ dplyr::case_when(primers == p_tag ~ p_notag, .default = .)
                    )
                )
            
            message(paste0("Primer parameter '",p_name,"' has been updated to '",p_notag,"' for 'pcr_primer' value '",p_tag,"'."))
        
        } else if ( stringr::str_detect(p_value, "^\\[\\S+\\][^\\[\\]\\;]+;\\[\\S+\\][^\\[\\]\\;]+(\\[\\S+\\][^\\[\\]\\;]+)*$") ) { # [A]1;B[2](;[C]3) -- two or more values, each with a primers tag
        
            # split value into parts using semi-colons
            p_value_split <-  stringr::str_split_1(p_value, ";")
            
            # extract tags from each part
            p_tags <- stringr::str_extract(p_value_split, "\\[(\\S+)\\]", group = 1)
            # extract true values from each part
            p_notags <- stringr::str_extract(p_value_split, "\\[\\S+\\](\\S+)", group = 1)
            
            # check all tags are valid primers values
            if ( !all(p_tags %in% primers_vec )) {
                stop(paste0(ppe, "The value of parameter '--pp_",p_name,"' does not have one or more valid primers tags: '", p_value,"'"))
            }
            
            # create tibble of parameter values ID'd by primers value
            names(p_notags) <- p_tags
            p_tibble <- tibble::enframe(p_notags, name = "primers", value = p_name)
            
            # update primer_params tibble with new values using primers as key
            primer_params_up <- 
                primer_params_up %>%
                dplyr::rows_update(
                .,
                p_tibble, 
                by = "primers",
                unmatched = "error"
                )
            
            message(paste0("Primer parameter '",p_name,"' has been updated to '",stringr::str_flatten(p_notags, " / "),"' for 'pcr_primer' values '",stringr::str_flatten(p_tags, " / "),"', respectively."))
        
        } else {
            stop(paste0(ppe, "The value of parameter '--pp_",p_name,"' is misformatted: '", p_value,"'"))
        }
    }
} else {
    message("No primer parameters replaced.")
}

# convert phmm, idtaxa_db and ref_fasta paths to absolute paths
primer_params_up <- 
    primer_params_up %>% 
    dplyr::mutate(
        dplyr::across(
            c(phmm, idtaxa_db, ref_fasta), 
            ~ dplyr::case_when(
                is.na(.) | . == "NA" ~ ., # keep as NA if NA
                stringr::str_starts(., "/") ~ ., # if already absolute path, leave it be
                stringr::str_starts(., "\\./") ~ stringr::str_replace(., "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(., "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    ) 

## parameter validation

# if `idtaxa_db` is not provided and `params.train_idtaxa` is not true, throw error
if ( any(primer_params_up$idtaxa_db %in% c("null",NA) ) && params.train_idtaxa %in% c("null","false") ) {
    stop (paste0(prpe, "One or more values of 'idtaxa_db' are missing and 'params.train_idtaxa' is set to 'null' or 'false'."))
}

# check locus, for_primer_seq and rev_primer_seq are not NA
if (any(is.na(primer_params_up$locus))) {
    stop(paste0(prpe, "One or more values of 'locus' are missing: '",stringr::str_flatten(primer_params_up$locus, " / "),"'."))
}
if (any(is.na(primer_params_up$for_primer_seq))) {
    stop(paste0(prpe, "One or more values of 'for_primer_seq' are missing: '",stringr::str_flatten(primer_params_up$for_primer_seq, " / "),"'."))
}
if (any(is.na(primer_params_up$rev_primer_seq))) {
    stop(paste0(prpe, "One or more values of 'rev_primer_seq' are missing: '",stringr::str_flatten(primer_params_up$rev_primer_seq, " / "),"'."))
}

# check for_primer_seq and rev_primer_seq are correctly formatted
if (!all(stringr::str_detect(primer_params_up$for_primer_seq, "^[GATCRYMKSWHBVDNIgatcrymkswhbvdni]+$"))){
    stop(paste0(prpe, "One or more values of 'for_primer_seq' contains non-IUPAC characters: '", stringr::str_flatten(primer_params_up$for_primer_seq, " / "),"'."))
}
if (!all(stringr::str_detect(primer_params_up$rev_primer_seq, "^[GATCRYMKSWHBVDNIgatcrymkswhbvdni]+$"))){
    stop(paste0(prpe, "One or more values of 'rev_primer_seq' contains non-IUPAC characters: '", stringr::str_flatten(primer_params_up$rev_primer_seq, " / "),"'."))
}

# check that still none of the primer_params values contain spaces or commas 
if (!all(stringr::str_detect(primer_params_up %>% unlist %>% unname, "^[^\\s,;]+$") %>% .[!is.na(.)] )){
    stop(paste0(prpe, "Primer parameters must not contain spaces, commas or semi-colons."))
}

# check max_primer_mismatch is an integer >= 0
if (any(primer_params_up$max_primer_mismatch %>% as.integer() %>% is.na()) || any(primer_params_up$max_primer_mismatch %>% as.integer() < 0)){
    stop(paste0(prpe, "'max_primer_mismatch' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$max_primer_mismatch, " / "), "'."))
}

# check read_min_length is an integer >= 0
if (any(primer_params_up$read_min_length %>% as.integer() %>% is.na()) || any(primer_params_up$read_min_length %>% as.integer() < 0)){
    stop(paste0(prpe, "'read_min_length' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$read_min_length, " / "), "'."))
}

# check read_max_length is an integer >= 1 or 'Inf'
for (i in 1:length(primer_params_up$read_max_length) ) {
    test_val <- primer_params_up$read_max_length[i]
    if ( test_val != "Inf" ){
        if ( test_val %>% as.integer %>% is.na || test_val %>% as.integer < 1 ){
            stop(paste0(prpe, "'read_max_length' must be an integer >=1, or the string 'Inf': '",stringr::str_flatten(primer_params_up$read_max_length, " / "), "'."))
        }
    }
}

# check read_max_ee is an integer >= 0
if (any(primer_params_up$read_max_ee %>% as.integer() %>% is.na()) || any(primer_params_up$read_max_ee %>% as.integer() < 0)){
    stop(paste0(prpe, "'read_max_ee' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$read_max_ee, " / "), "'."))
}

# check read_trunc_length is an integer >= 1
if (any(primer_params_up$read_trunc_length %>% as.integer() %>% is.na()) || any(primer_params_up$read_trunc_length %>% as.integer() < 0)){
    stop(paste0(prpe, "'read_trunc_length' must be an integer >= 1: '",stringr::str_flatten(primer_params_up$read_trunc_length, " / "), "'."))
}

# check read_trim_left is an integer >= 0
if (any(primer_params_up$read_trim_left %>% as.integer() %>% is.na()) || any(primer_params_up$read_trim_left %>% as.integer() < 0)){
    stop(paste0(prpe, "'read_trim_left' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$read_trim_left, " / "), "'."))
}

# check read_trim_right is an integer >= 0
if (any(primer_params_up$read_trim_right %>% as.integer() %>% is.na()) || any(primer_params_up$read_trim_right %>% as.integer() < 0)){
    stop(paste0(prpe, "'read_trim_right' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$read_trim_right, " / "), "'."))
}

# check asv_min_length is an integer >= 1
if (any(primer_params_up$asv_min_length %>% as.integer() %>% is.na()) || any(primer_params_up$asv_min_length %>% as.integer() < 0)){
    stop(paste0(prpe, "'asv_min_length' must be an integer >= 1: '",stringr::str_flatten(primer_params_up$asv_min_length, " / "), "'."))
}

# check asv_max_length is an integer >= 1
if (any(primer_params_up$asv_max_length %>% as.integer() %>% is.na()) || any(primer_params_up$asv_max_length %>% as.integer() < 0)){
    stop(paste0(prpe, "'asv_max_length' must be an integer >= 1: '",stringr::str_flatten(primer_params_up$asv_max_length, " / "), "'."))
}

# check concat_unmerged is a logical string
if ( !all(stringr::str_detect(toupper(primer_params_up$concat_unmerged), "^(TRUE|T|FALSE|F)$")) ){
    stop(paste0(prpe, "'concat_unmerged' must be one of 'TRUE', 'T', 'FALSE' or 'F' (or lowercase versions): '",stringr::str_flatten(primer_params_up$concat_unmerged, " / "), "'."))
}

# check genetic code is a valid code from Biostrings::GENETIC_CODE_TABLE, or NA
gc_valid <- c(Biostrings::GENETIC_CODE_TABLE$name2, Biostrings::GENETIC_CODE_TABLE$id, NA) %>% unique()
if ( !all(primer_params_up$genetic_code %in% gc_valid) ){
    stop(paste0(prpe, "'genetic_code' must be a 'name2' or 'id' value from Biostrings::GENETIC_CODE_TABLE: '",stringr::str_flatten(primer_params_up$genetic_code, " / "),"'."))
}

# check coding is a logical string
if ( !all(stringr::str_detect(toupper(primer_params_up$coding), "^(TRUE|T|FALSE|F)$")) ){
    stop(paste0(prpe, "'coding' must be one of 'TRUE', 'T', 'FALSE' or 'F' (or lowercase versions): '",stringr::str_flatten(primer_params_up$coding, " / "), "'."))
}

## check phmm, idtaxa_db and ref_fasta are readable files
# phmm path check 
check_paths <- primer_params_up$phmm %>% unlist() # check phmm paths
for(i in seq_along(check_paths)){ if (!is.na(check_paths[i])) {
    message (paste0("Checking 'phmm' path '",check_paths[i],"'..."))
    assertthat::is.readable(check_paths[i])} 
}

# idtaxa_db path check
check_paths <- primer_params_up$idtaxa_db %>% unlist() # check idtaxa_db paths
for(i in seq_along(check_paths)){ if (!is.na(check_paths[i])) {
    message (paste0("Checking 'idtaxa_db' path '",check_paths[i],"'..."))
    assertthat::is.readable(check_paths[i])} 
}

# ref_fasta path check
check_paths <- primer_params_up$ref_fasta %>% unlist() # check ref_fasta paths
for(i in seq_along(check_paths)){ 
    message (paste0("Checking 'ref_fasta' path '",check_paths[i],"'..."))
    assertthat::is.readable(check_paths[i]) 
}

##
# check idtaxa_confidence is an integer 0-100
if (any(primer_params_up$idtaxa_confidence %>% as.integer() %>% is.na()) || !all(primer_params_up$idtaxa_confidence %>% as.integer() %>% between(., 0, 100))){
    stop(paste0(prpe, "'idtaxa_confidence' must be an integer between 0-100: '",stringr::str_flatten(primer_params_up$idtaxa_confidence, " / "), "'."))
}

# check run_blast is a logical string
if ( !all(stringr::str_detect(toupper(primer_params_up$run_blast), "^(TRUE|T|FALSE|F)$")) ){
    stop(paste0(prpe, "'run_blast' must be one of 'TRUE', 'T', 'FALSE' or 'F' (or lowercase versions): '",stringr::str_flatten(primer_params_up$run_blast, " / "), "'."))
}

# check blast_min_identity is an integer 0-100
if (any(primer_params_up$blast_min_identity %>% as.integer() %>% is.na()) || !all(primer_params_up$blast_min_identity %>% as.integer() %>% between(., 0, 100))){
    stop(paste0(prpe, "'blast_min_identity' must be an integer between 0-100: '",stringr::str_flatten(primer_params_up$blast_min_identity, " / "), "'."))
}

# check blast_min_coverage is an integer 0-100
if (any(primer_params_up$blast_min_coverage %>% as.integer() %>% is.na()) || !all(primer_params_up$blast_min_coverage %>% as.integer() %>% between(., 0, 100))){
    stop(paste0(prpe, "'blast_min_coverage' must be an integer between 0-100: '",stringr::str_flatten(primer_params_up$blast_min_coverage, " / "), "'."))
}

# NOTE: target_* fields should already be validated as lacking spaces or commas -- invalid taxa names should be ignored at this stage

# check min_sample_reads is an integer >= 0
if (any(primer_params_up$min_sample_reads %>% as.integer() %>% is.na()) || any(primer_params_up$min_sample_reads %>% as.integer() < 0)){
    stop(paste0(prpe, "'min_sample_reads' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$min_sample_reads, " / "), "'."))
}

# check min_taxa_reads is an integer >= 0
if (any(primer_params_up$min_taxa_reads %>% as.integer() %>% is.na()) || any(primer_params_up$min_taxa_reads %>% as.integer() < 0)){
    stop(paste0(prpe, "'min_taxa_reads' must be an integer >= 0: '",stringr::str_flatten(primer_params_up$min_taxa_reads, " / "), "'."))
}

# check min_taxa_ra is a number between 0-1
if (any(primer_params_up$min_taxa_ra %>% as.numeric() %>% is.na()) || !all(primer_params_up$min_taxa_ra %>% as.numeric() %>% between(., 0, 1)) ){
    stop(paste0(prpe, "'min_taxa_ra' must be a number between 0-1: '",stringr::str_flatten(primer_params_up$min_taxa_ra, " / "), "'."))
}

## write parsed primer_params to file
readr::write_csv(primer_params_up, "primer_params_parsed.csv")

# stop("stopped manually")