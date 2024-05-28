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
    "projectDir",
    "samplesheet",
    "loci_params"
)
lapply(nf_vars, nf_var_check)

### run code

samplesheet_df <- readr::read_csv(paste0(projectDir,"/",samplesheet), show_col_types = F) # read in samplesheet

loci_params_df <- readr::read_csv(paste0(projectDir,"/",loci_params), show_col_types = F) # read in loci parameters


### validation of samplesheet content

for (i in 1:length(samplesheet_df$sample_id)) { # check that sample_name is a subset of sample_id
    if (!stringr::str_detect(samplesheet_df$sample_id[i], samplesheet_df$sample_name[i])) {
        stop(paste0("INPUT ERROR: 'sample_name' (",samplesheet_df$sample_name[i],") is not a subset of 'sample_id' (", samplesheet_df$sample_id[i],")!\n\t'sample_id' must be a unique combination of 'fcid' and 'sample_name'."))
    }
}

for (i in 1:length(samplesheet_df$sample_id)) { # check that sample_id starts with fcid
    if (!stringr::str_starts(samplesheet_df$sample_id[i], samplesheet_df$fcid[i])) {
        stop(paste0("INPUT ERROR: 'sample_id' (", samplesheet_df$sample_id[i],") does not begin with 'fcid' (",samplesheet_df$fcid[i],")!\n\t'sample_id' must be a unique combination of 'fcid' and 'sample_name'."))
    }
}

if (any(duplicated(samplesheet_df$sample_id))) { # check that sample_id contains only unique values
    stop("INPUT ERROR: All sample IDs in samplesheet must be unique!") 
    } 

# check that all primer names and sequences are identical for each sample
if (!any(duplicated(samplesheet_df$pcr_primers))) { stop("INPUT ERROR: 'pcr_primer' varies between samples!\n\t'pcr_primer' contents must be the same for all samples.") } 
if (!any(duplicated(samplesheet_df$for_primer_seq))) { stop("INPUT ERROR: 'for_primer_seq' varies between samples!\n\t'for_primer_seq' contents must be the same for all samples.") } 
if (!any(duplicated(samplesheet_df$rev_primer_seq))) { stop("INPUT ERROR: 'rev_primer_seq' varies between samples!\n\t'rev_primer_seq' contents must be the same for all samples.") } 

# check that same number of target genes and primer names are defined as primer sequences
if (
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) |
    !identical(stringr::str_count(samplesheet_df$target_gene, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("INPUT ERROR: Mismatch between number of supplied target genes and primer sequences!") }
if (
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$for_primer_seq, ";")) |
    !identical(stringr::str_count(samplesheet_df$pcr_primers, ";"), stringr::str_count(samplesheet_df$rev_primer_seq, ";"))
    ) { stop("INPUT ERROR: Mismatch between number of supplied primer pair names and primer sequences!") }

### parse and/or detect read paths

if ( "read_dir" %in% colnames(samplesheet_df) ) { # if "read_dir" column exists in samplesheet
    if ("fwd" %in% colnames(samplesheet_df) | "rev" %in% colnames(samplesheet_df) ) { # if "fwd" or "rev" column exists as well as "read_dir", throw error
        stop ("INPUT ERROR: 'read_dir' and one or both of 'fwd' and 'rev' should not be specified together.\n\tSpecify read file locations using only one of: the read directory path ('read_dir'), or direct read paths ('fwd' and 'rev').")
    } else { # try to find reads
        # convert read directory to absolute path
        samplesheet_df <- samplesheet_df %>% 
            dplyr::mutate(
                read_dir = dplyr::case_when(
                    stringr::str_starts(read_dir, "/") ~ read_dir, # if already absolute path, leave it be
                    stringr::str_starts(read_dir, "\\./") ~ stringr::str_replace(read_dir, "^\\.", projectDir),  # replace "." with projectDir to produce absolute path
                    .default = stringr::str_replace(read_dir, "^", paste0(projectDir,"/")) # else lead with projectDir to produce absolute path
                )
            )
        reads_list = list() # create empty list
        for (i in 1:length(samplesheet_df$read_dir)) { # loop through rows of samplesheet
            i_readfiles <- list.files( # find full paths of the reads matching sample_id
                path = samplesheet_df$read_dir[i],
                pattern = samplesheet_df$sample_id[i],
                full.names = T, 
                recursive = T
                ) %>% unlist()
            # check exactly two read files are found
            if ( length(i_readfiles) != 2 ) { 
                stop (paste0("ERROR: Found ",length(i_readfiles)," read files matching '",samplesheet_df$sample_id[i],"' in '",samplesheet_df$read_dir[i],"' when 2 were expected.\n\tCheck you have filled out samplesheet correctly."))  
            }
            if (exists("params.extension")) { # using params.extension if supplied
                message(paste0("Using 'params.extension (",params.extension,") to find reads in supplied reads directories ('read_dir')."))
                # check read files match params.extension
                if (any(stringr::str_detect(i_readfiles, pattern = paste0(params.extension,"$")))) { stop("Error: Found read files do not match provided extension--check samplesheet.")}
                # sort read files by text after extension
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
    } 
} else { # if "read_dir" doesn't exist...
    if ("fwd" %in% colnames(samplesheet_df) & "rev" %in% colnames(samplesheet_df) ) { # if 'fwd' and 'rev' columns both exist...
        # convert read paths to absolute paths
        samplesheet_df <- samplesheet_df %>% 
            dplyr::mutate(
                fwd = dplyr::case_when(
                    stringr::str_starts(fwd, "/") ~ fwd, # if already absolute path, leave it be
                    stringr::str_starts(fwd, "\\./") ~ stringr::str_replace(fwd, "^\\.", projectDir),  # replace "." with projectDir to produce absolute path
                    .default = stringr::str_replace(fwd, "^", paste0(projectDir,"/")) # else append projectDir to front to produce absolute path
                ),
                rev = dplyr::case_when(
                    stringr::str_starts(rev, "/") ~ rev, # if already absolute path, leave it be
                    stringr::str_starts(rev, "\\./") ~ stringr::str_replace(rev, "^\\.", projectDir),  # remove "." and lead with projectDir to produce absolute path
                    .default = stringr::str_replace(rev, "^", paste0(projectDir,"/")) # else lead with projectDir to produce absolute path
                )
            )
        
    } else {
        stop ("INPUT ERROR: one of 'fwd' or 'rev' is not present, with 'read_dir' not present!")
    }
}

# glimpse(samplesheet_df)

## check that read file paths are readable
check_paths <- samplesheet_df$fwd %>% unlist() # check fwd paths
for(i in seq_along(check_paths)){ assertthat::is.readable(check_paths[i]) }
check_paths <- samplesheet_df$rev %>% unlist() # check rev paths
for(i in seq_along(check_paths)){ assertthat::is.readable(check_paths[i]) }

# check that read file paths are unique
### TODO: make this a more informative error -- say which files are duplicated
if (!any(duplicated(samplesheet_df$fwd))) {stop ("INPUT ERROR: At least two samples share the same forward read file in the samplesheet!")}
if (any(duplicated(samplesheet_df$rev))) {stop ("INPUT ERROR: At least two samples share the same reverse read files in the samplesheet!")}

## write new samplesheet to file


### validate loci_params content

# add missing columns with default values

# check pcr_primers column contains unique values

# check phmm, idtaxa_db and ref_fasta are readable files






### split samplesheets into loci and join to loci_params


stop(" *** stopped manually *** ") ##########################################
