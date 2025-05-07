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
    "loci_params",
    "seq_type",
    "paired"
)
lapply(nf_vars, nf_var_check)

### process variables 

# ## make relative paths absolute
# # for samplesheet path
# if ( stringr::str_starts(samplesheet, "\\./") ) {
#     samplesheet <- stringr::str_replace(samplesheet, "\\.", launchDir)
# } else if ( !stringr::str_starts(samplesheet, "/") ) {
#     samplesheet <- paste0(launchDir,"/",samplesheet)
# } else {
#     samplesheet <- samplesheet
# }

# # for loci_params path
# if ( stringr::str_starts(loci_params, "\\./") ) {
#     loci_params <- stringr::str_replace(loci_params, "\\.", launchDir)
# } else if ( !stringr::str_starts(loci_params, "/") ) {
#     loci_params <- paste0(launchDir,"/",loci_params)
# } else {
#     loci_params <- loci_params
# }

### run code

# read in samplesheet
samplesheet_df <- readr::read_csv(samplesheet, show_col_types = F, skip_empty_rows = TRUE, na = c("", "NA", " ")) %>%
    dplyr::select(!tidyselect::starts_with("...")) %>% # remove any blank columns that have been named "...X"
    dplyr::filter(dplyr::if_any(tidyselect::everything(), ~ !is.na(.))) # remove completely blank lines

# pseudorandomly subsample samplesheet if params.subsample is defined
# this is done per pcr_primer x fcid combination so expected combinations are (likely) retained
set.seed(1)
if (params.subsample != "null") {
    samplesheet_df <- samplesheet_df %>% 
    dplyr::group_by(pcr_primers, fcid) %>% 
    dplyr::slice_sample(n = as.numeric(params.subsample), replace = F) %>%
    dplyr::ungroup()
}

# read in loci parameters
loci_params_df <- readr::read_csv(loci_params, show_col_types = F, skip_empty_rows = TRUE, na = c("", "NA", " ")) %>%
    dplyr::select(!tidyselect::starts_with("...")) %>% # remove any blank columns that have been named "...X"
    dplyr::filter(dplyr::if_any(tidyselect::everything(), ~ !is.na(.))) # remove completely blank lines

### validation of samplesheet content

## add missing columns if not present, with default values
# make default samplesheet
samplesheet_default <- tibble::tibble( # create columns with NA values in a single row
    sample_id = NA_character_,
    sample_name = NA_character_,
    extraction_rep = NA_character_,
    amp_rep = NA_character_,
    client_name = NA_character_,
    experiment_name = NA_character_,
    sample_type = NA_character_,
    collection_method = NA_character_,
    collection_location = NA_character_,
    latitude = NA_integer_,
    longitude = NA_integer_,
    environment = NA_character_,
    collection_date = NA_character_,
    operator_name = NA_character_,
    description = NA_character_,
    assay = NA_character_,
    extraction_method = NA_character_,
    amp_method = NA_character_,
    target_gene = NA_character_,
    pcr_primers = NA_character_,
    for_primer_seq = NA_character_,
    rev_primer_seq = NA_character_,
    index_plate = NA_character_,
    index_well = NA_character_,
    i7_index_id = NA_character_,
    i7_index = NA_character_,
    i5_index_id = NA_character_,
    i5_index = NA_character_,
    seq_platform = NA_character_,
    fcid = NA_character_,
    for_read_length = NA_integer_,
    rev_read_length = NA_integer_,
    seq_run_id = NA_character_,
    seq_id = NA_character_,
    seq_date = NA_character_,
    analysis_method = NA_character_,
    notes = NA_character_
) 

# combine defaults with input samplesheet
samplesheet_df <- new_bind(samplesheet_default %>% filter(FALSE), samplesheet_df)

# stop(" *** stopped manually *** ") ##########################################

## validate contents of columns
for (i in 1:length(samplesheet_df$sample_id)) { # check that sample_name is a subset of sample_id
    if (!stringr::str_detect(samplesheet_df$sample_id[i], samplesheet_df$sample_name[i])) {
        stop(paste0("SAMPLESHEET ERROR: 'sample_name' (",samplesheet_df$sample_name[i],") is not a subset of 'sample_id' (", samplesheet_df$sample_id[i],")!\n\t'sample_id' must be a unique combination of 'fcid' and 'sample_name'."))
    }
}

for (i in 1:length(samplesheet_df$sample_id)) { # check that sample_id starts with fcid
    if (!stringr::str_starts(samplesheet_df$sample_id[i], samplesheet_df$fcid[i])) {
        stop(paste0("SAMPLESHEET ERROR: 'sample_id' (", samplesheet_df$sample_id[i],") does not begin with 'fcid' (",samplesheet_df$fcid[i],")!\n\t'sample_id' must be a unique combination of 'fcid' and 'sample_name'."))
    }
}

if (any(duplicated(samplesheet_df$sample_id))) { # check that sample_id contains only unique values
    stop("SAMPLESHEET ERROR: Some sample IDs are repeated!\n\tAll 'sample_id' values in samplesheet must be unique.") 
    } 

# check that all primer names and sequences are identical for each sample
if (nrow(samplesheet_df) != 1) { # checks fail if samplesheet only has one row
    if (!any(duplicated(samplesheet_df$pcr_primers))) { stop("SAMPLESHEET ERROR: 'pcr_primers' varies between samples!\n\t'pcr_primers' contents must be the same for all samples.") } 
    if (!any(duplicated(samplesheet_df$for_primer_seq))) { stop("SAMPLESHEET ERROR: 'for_primer_seq' varies between samples!\n\t'for_primer_seq' contents must be the same for all samples.") } 
    if (!any(duplicated(samplesheet_df$rev_primer_seq))) { stop("SAMPLESHEET ERROR: 'rev_primer_seq' varies between samples!\n\t'rev_primer_seq' contents must be the same for all samples.") } 
}
# check that same number of target genes and primer names are defined as primer sequences
if (
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) |
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("SAMPLESHEET ERROR: Mismatch between number of supplied target genes and primer sequences!") }
if (
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) |
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("SAMPLESHEET ERROR: Mismatch between number of supplied primer pair names and primer sequences!") }

# check if disallowed read column combinations are detected
if ( paired == "true" & "single" %in% colnames(samplesheet_df) ) {
    stop( "SAMPLESHEET ERROR: If data is paired-end, 'single' field cannot be used in the samplesheet.\n\tUse either 'fwd' and 'rev', or 'read_dir' to specify data locations." )
} 

if ( paired == "false" & ( "fwd" %in% colnames(samplesheet_df) | "rev" %in% colnames(samplesheet_df) ) ) {
    stop( "SAMPLESHEET ERROR: If data is single-end, 'fwd' or 'rev' fields cannot be used in the samplesheet.\n\tUse either 'single', or 'read_dir' to specify data locations." )
} 

if ( 
    ( "fwd" %in% colnames(samplesheet_df) & "single" %in% colnames(samplesheet_df) ) | 
    ( "rev" %in% colnames(samplesheet_df) & "single" %in% colnames(samplesheet_df) ) | 
    ( "read_dir" %in% colnames(samplesheet_df) & "fwd" %in% colnames(samplesheet_df) ) | 
    ( "read_dir" %in% colnames(samplesheet_df) & "rev" %in% colnames(samplesheet_df) ) | 
    ( "read_dir" %in% colnames(samplesheet_df) & "single" %in% colnames(samplesheet_df) ) 
    ) { 
    stop ("SAMPLESHEET ERROR: Disallowed combination of read file columns.\n\tSpecify read file locations using only one of: the read directory path ('read_dir'), paired-end direct read paths ('fwd' and 'rev'), or single-end direct read paths ('single').") 
    }

### parse and/or detect read paths
if ( "read_dir" %in% colnames(samplesheet_df) & paired == "true" ) { # if "read_dir" column exists in samplesheet
    ## try to find reads
    # convert read directory to absolute path
    samplesheet_df <- samplesheet_df %>% 
        dplyr::mutate(
            read_dir = dplyr::case_when(
                stringr::str_starts(read_dir, "/") ~ read_dir, # if already absolute path, leave it be
                stringr::str_starts(read_dir, "\\./") ~ stringr::str_replace(read_dir, "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(read_dir, "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    reads_list = list() # create empty list
    for (i in 1:length(samplesheet_df$read_dir)) { # loop through rows of samplesheet
        i_readfiles <- list.files( # find full paths of files matching sample_id with a fastQ extension
            path = samplesheet_df$read_dir[i],
            pattern = paste0(samplesheet_df$sample_id[i],"(?:_[^\\s/]+)?\\.f(ast)?q(\\.gz)?$"),
            full.names = T, 
            recursive = T
            ) %>% unlist()
        # check exactly two read files are found
        if ( length(i_readfiles) != 2 ) { 
            stop (paste0("SAMPLESHEET ERROR: Found ",length(i_readfiles)," read files matching '",samplesheet_df$sample_id[i],"' in '",samplesheet_df$read_dir[i],"' when 2 were expected.\n\tCheck you have filled out samplesheet correctly."))  
        }
        if (params.extension != "null") { # using params.extension if supplied
            message(paste0("Using 'params.extension' (",params.extension,") to find read files in supplied directories ('read_dir')."))
            # check read files match params.extension
            if (any(stringr::str_detect(i_readfiles, pattern = paste0(params.extension,"$")))) { stop("SAMPLESHEET ERROR: Read files found in 'read_dir' do not match 'params.extension' file extension--check samplesheet and pipeline parameters.")}
            # sort read files by text after extension
            ### TODO: do this
        } else { # not using params.extension
            # sort read files by natural order (ie. 1 before 2; F before R)
            i_readfiles <- stringr::str_sort(i_readfiles)
        }
        
        reads_list[[i]] <- c(samplesheet_df$sample_id[i], i_readfiles) # append to list
    }
    reads_df <- do.call(rbind, reads_list) # combine list into dataframe
    colnames(reads_df) <- c("sample_id","fwd","rev") # name columns
    reads_df <- reads_df %>% tibble::as_tibble() # convert to tibble
    # add read paths to samplesheet
    samplesheet_df <- samplesheet_df %>% 
        dplyr::left_join(., reads_df, by = "sample_id") 

} else if ( "fwd" %in% colnames(samplesheet_df) & "rev" %in% colnames(samplesheet_df) & paired == "true") { # if 'fwd' and 'rev' columns both exist...
    # convert read paths to absolute paths
    samplesheet_df <- samplesheet_df %>% 
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
} else if ( "read_dir" %in% colnames(samplesheet_df) & paired == "false" ) { # if "read_dir" column exists in samplesheet
    ## try to find reads
    # convert read directory to absolute path
    samplesheet_df <- samplesheet_df %>% 
        dplyr::mutate(
            read_dir = dplyr::case_when(
                stringr::str_starts(read_dir, "/") ~ read_dir, # if already absolute path, leave it be
                stringr::str_starts(read_dir, "\\./") ~ stringr::str_replace(read_dir, "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(read_dir, "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    reads_list = list() # create empty list
    for (i in 1:length(samplesheet_df$read_dir)) { # loop through rows of samplesheet
        i_readfiles <- list.files( # find full paths of files matching sample_id with a fastQ extension
            path = samplesheet_df$read_dir[i],
            pattern = paste0(samplesheet_df$sample_id[i],"(?:_[^\\s/]+)?\\.f(ast)?q(\\.gz)?$"),
            full.names = T, 
            recursive = T
            ) %>% unlist()
        # check exactly one read file is found per sample
        if ( length(i_readfiles) != 1 ) { 
            stop (paste0("SAMPLESHEET ERROR: Found ",length(i_readfiles)," read files matching '",samplesheet_df$sample_id[i],"' in '",samplesheet_df$read_dir[i],"' when 1 was expected.\n\tCheck you have filled out samplesheet correctly."))  
        }
        if (params.extension != "null") { # using params.extension if supplied
            message(paste0("Using 'params.extension' (",params.extension,") to find read files in supplied directories ('read_dir')."))
            # check read files match params.extension
            if (any(stringr::str_detect(i_readfiles, pattern = paste0(params.extension,"$")))) { stop("SAMPLESHEET ERROR: Read files found in 'read_dir' do not match 'params.extension' file extension--check samplesheet and pipeline parameters.")}
            # sort read files by text after extension
            ### TODO: do this
        } else { # not using params.extension
            # sort read files by natural order (ie. 1 before 2; F before R)
            i_readfiles <- stringr::str_sort(i_readfiles)
        }
        
        reads_list[[i]] <- c(samplesheet_df$sample_id[i], i_readfiles) # append to list
    }
    reads_df <- do.call(rbind, reads_list) # combine list into dataframe
    colnames(reads_df) <- c("sample_id","single") # name columns
    reads_df <- reads_df %>% tibble::as_tibble() # convert to tibble
    # add read paths to samplesheet
    samplesheet_df <- samplesheet_df %>% 
        dplyr::left_join(., reads_df, by = "sample_id") 

} else if (( "single" %in% colnames(samplesheet_df) & paired == "false")) { # if single reads only
# convert read paths to absolute paths
    samplesheet_df <- samplesheet_df %>% 
        dplyr::mutate(
            across(
                c(single),
                ~ dplyr::case_when(
                    stringr::str_starts(., "/") ~ ., # if already absolute path, leave it be
                    stringr::str_starts(., "\\./") ~ stringr::str_replace(., "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                    .default = stringr::str_replace(., "^", paste0(launchDir,"/")) # else append launchDir to front to produce absolute path
                )
            )
        )
} else {
    stop ("SAMPLESHEET ERROR: Disallowed combination of samplesheet read columns and 'paired' pipeline parameter.")
}



# check that read file paths are unique
### TODO: make this a more informative error -- say which files are duplicated
if (paired == "true") {
    if (any(duplicated(samplesheet_df$fwd))) {stop ("SAMPLESHEET ERROR: At least two samples share the same forward read file in the samplesheet!")}
    if (any(duplicated(samplesheet_df$rev))) {stop ("SAMPLESHEET ERROR: At least two samples share the same reverse read file in the samplesheet!")}
} else if ( paired == "false" ) {
    if (any(duplicated(samplesheet_df$single))) {stop ("SAMPLESHEET ERROR: At least two samples share the same read file in the samplesheet!")}
} else {
    stop(paste0("PIPELINE PARAMETER ERROR: 'paired' = '", paired, "' but must be 'true' or 'false'."))
}


## write parsed samplesheet to file
readr::write_csv(samplesheet_df, "samplesheet_parsed.csv")


### validate loci_params content

# add missing columns with default values
default_params <- tibble::tibble(
    pcr_primers = NA_character_,
    target_gene = NA_character_,
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

loci_params_df <- new_bind(default_params %>% filter(FALSE), loci_params_df) %>%
  dplyr::mutate(across(everything(), as.character))%>%
  dplyr::select(colnames(default_params))

### replace parameters from '--loci_params' .csv file with those from the '--lp_*' flags
# vector of loci parameters
lp_vec <- c(
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
)

# loop through parameters, changing those with lp_* not set to "null"
for ( i in 1:length(lp_vec) ) {
    lp <- lp_vec[i]
    p <- get(paste0("params.lp_",lp))
    if ( p == "null" ) { NULL } else {
        if ( p %>% stringr::str_detect(pattern = "^<.+>.+$")) {
        ## if multi-locus values given
        # get the loci names
        p_names <- stringr::str_extract_all(p, pattern = "(?<=\\<).+?(?=\\>)") %>% unlist
        # get the new parameter values
        p_vals <- stringr::str_extract_all(p, pattern = "(?<=\\>).+?((?=\\<)|$)") %>% unlist
        # create a tibble of the new parameter values
        p_df <- data.frame(p_names, p_vals) %>% 
            tibble::as_tibble() %>% 
            dplyr::mutate(across(everything(), as.character)) %>%
            dplyr::rename("pcr_primers" = p_names, !!lp := p_vals) 
        # update the original df with the new parameter values per locus
        loci_params_df <- dplyr::rows_update(loci_params_df, p_df, by = "pcr_primers")
        } else {
            ## single value given, replace all old values with new value
            loci_params_df <- loci_params_df %>%
                dplyr::mutate_at( vars(lp), function(x) (x = p) )
        }
    }
}

# if `idtaxa_db` is not provided and `params.train_idtaxa` is not true, throw error
if ( any(loci_params_df$idtaxa_db %in% c("null",NA) ) && params.train_idtaxa == "null" ) {
    stop ("LOCI_PARAMS ERROR: One or more values of 'idtaxa_db' are missing and 'params.train_idtaxa' is set to 'null'.")
    }

# check pcr_primers column contains only unique values
if (any(duplicated(loci_params_df$pcr_primers))) {
    stop ("LOCI_PARAMS ERROR: 'pcr_primers' values are not unique!")
    }

# create split samplesheet for checking against loci_params
samplesheet_split_check <- samplesheet_df %>% 
    tidyr::separate_longer_delim(c(pcr_primers, target_gene, for_primer_seq, rev_primer_seq), delim = ";") 

# check pcr_primers values match those in samplesheet
if (!setequal(samplesheet_split_check$pcr_primers %>% unique(), loci_params_df$pcr_primers)) {
    stop("LOCI_PARAMS ERROR: Samplesheet and loci_params do not share the same values of 'pcr_primers'!")
}

# check target_gene values match those in samplesheet
if (!setequal(samplesheet_split_check$target_gene %>% unique(), loci_params_df$target_gene)) {
    stop("LOCI_PARAMS ERROR: Samplesheet and loci_params do not share the same values of 'target_gene'!")
}

# convert run_blast, coding and concat_unmerged to logical values
loci_params_df <- loci_params_df %>%
    dplyr::mutate( ### NOTE: this assumes all columns are present in the data frame, but this might not always be true
        run_blast = as.logical(run_blast),
        coding = as.logical(coding),
        concat_unmerged = dplyr::case_when(
            concat_unmerged == "NA" ~ FALSE,
            concat_unmerged %in% c("TRUE", "true") ~ TRUE,
            concat_unmerged %in% c("FALSE", "false") ~ FALSE,
            .default = FALSE
        )
    )

# convert phmm, idtaxa_db and ref_fasta paths to absolute paths
loci_params_df <- loci_params_df %>% 
    dplyr::mutate(
        across(
            c(phmm, idtaxa_db, ref_fasta), 
            ~ dplyr::case_when(
                is.na(.) ~ ., # keep as NA if NA
                stringr::str_starts(., "/") ~ ., # if already absolute path, leave it be
                stringr::str_starts(., "\\./") ~ stringr::str_replace(., "^\\.", launchDir),  # replace "." with launchDir to produce absolute path
                .default = stringr::str_replace(., "^", paste0(launchDir,"/")) # else lead with launchDir to produce absolute path
            )
        )
    )

## check phmm, idtaxa_db and ref_fasta are readable files
# phmm path check 
check_paths <- loci_params_df$phmm %>% unlist() # check phmm paths
for(i in seq_along(check_paths)){ if (!is.na(check_paths[i])) {
    message (paste0("Checking 'phmm' path '",check_paths[i],"'..."))
    assertthat::is.readable(check_paths[i])} 
    }

# idtaxa_db path check
check_paths <- loci_params_df$idtaxa_db %>% unlist() # check idtaxa_db paths
for(i in seq_along(check_paths)){ if (!is.na(check_paths[i])) {
    message (paste0("Checking 'idtaxa_db' path '",check_paths[i],"'..."))
    assertthat::is.readable(check_paths[i])} 
    }

# ref_fasta path check
check_paths <- loci_params_df$ref_fasta %>% unlist() # check ref_fasta paths
for(i in seq_along(check_paths)){ assertthat::is.readable(check_paths[i]) }

## write parsed loci_parms to file
readr::write_csv(loci_params_df, "loci_params_parsed.csv")


### split samplesheet by primer and join to loci_params
samplesheet_loci_params <- samplesheet_df %>%  # split samdf loci-relevant columns across new rows, then join samdf and params  
    tidyr::separate_longer_delim(c(pcr_primers, for_primer_seq, rev_primer_seq, target_gene), delim = ";") %>% 
    dplyr::left_join(., loci_params_df, by = c("pcr_primers", "target_gene"))

readr::write_csv(samplesheet_loci_params, "samplesheet_loci_params.csv")

### split joined samplesheet into one per primer pair
split_slp <- split(samplesheet_loci_params, samplesheet_loci_params$pcr_primers) # split dfs by pcr_primers

for ( I in 1:length(split_slp)) { # assign new dfs to new variables
    new_df_name <- paste0(unique(split_slp[[I]]$pcr_primers),"__samplesheet")
    assign(
        paste0(unique(split_slp[[I]]$pcr_primers),"__samplesheet"),
        split_slp[[I]]
        )
    readr::write_csv( # print dfs inside work dir; maybe publish?
        x = get(new_df_name), 
        file = sprintf("%s.csv",new_df_name)
        )
}

# stop(" *** stopped manually *** ") ##########################################
