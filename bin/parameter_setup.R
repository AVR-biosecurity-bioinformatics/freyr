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


# Find all directories within data folder
if (!exists("params.data_folder")) { # if data_loc not defined, use "data"
    data_loc="data"
} else {
    data_loc = params.data_folder
}

data_loc_abs <- paste0(projectDir,"/",data_loc) # define absolute path for data directory
runs <- dir(data_loc_abs) # find subdirectories in data directory
SampleSheet <- list.files(paste0(data_loc_abs,"/", runs), pattern= "SampleSheet", full.names = TRUE) # list all sampleshets
runParameters <- list.files(paste0(data_loc_abs,"/", runs), pattern= "[Rr]unParameters.xml", full.names = TRUE) # list all run parameter files

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4") %>%
distinct()

# Check that sample_ids contain fcid, if not; attach
samdf <- samdf %>%
    dplyr::mutate(sample_id = case_when(
        !stringr::str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
        TRUE ~ sample_id
        )
    )

# Check that samples match samplesheet
fastqFs <- 
    purrr::map(list.dirs(data_loc_abs, recursive=FALSE),
                list.files, pattern="_R1_", full.names = TRUE) %>%
    unlist() %>%
    stringr::str_remove(pattern = "^(.*)\\/") %>%
    stringr::str_remove(pattern = "(?:.(?!_S))+$")

# Filter undetermined reads from sample sheet
fastqFs <- fastqFs[!stringr::str_detect(fastqFs, "Undetermined")]

# Check for fastq files that are missing from samplesheet
if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {
    warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") 
    }

# Check for sample_ids that dont have a corresponding fastq file
if (length(setdiff(samdf$sample_id, fastqFs)) > 0) {
    warning(paste0("The fastq file: ",
                    setdiff(samdf$sample_id, fastqFs),
                    " is missing, dropping from samplesheet \n")) 
    samdf <- samdf %>%
        filter(!sample_id %in% setdiff(samdf$sample_id, fastqFs))
}

### Add .fastq paths to the sample sheet in additional fwd and rev columns
# list all .fastq.gz files in data directory 
fastq_paths <- 
    list.files(data_loc_abs, pattern = "fastq.gz$|fq.gz$", recursive = T, full.names = T)

fastq_paths.fwd <- fastq_paths %>% stringr::str_subset(pattern = "_R1_") # list of fwd read paths
fastq_paths.rev <- fastq_paths %>% stringr::str_subset(pattern = "_R2_") # list of rev read paths
fastq_paths.base <- fastq_paths.fwd %>% stringr::str_extract(pattern = "([^\\/]+$)") # list of fwd read names 
fastq_paths.id <- fastq_paths.base %>% stringr::str_remove(pattern = "(?:.(?!_S))+$") # list of sample names from read names (everything before "_S")
fastq_paths.fcid <- fastq_paths.id %>% stringr::str_extract(pattern = "^([^_]+)") # fcid from sample name

fastq_paths.df <- # create data frame of read path info
    data.frame(fastq_paths.fwd, fastq_paths.rev, fastq_paths.base, fastq_paths.id, fastq_paths.fcid) %>%
    dplyr::rename(
        fwd = fastq_paths.fwd, 
        rev = fastq_paths.rev,
        base = fastq_paths.base,
        sample_id = fastq_paths.id,
        fcid = fastq_paths.fcid
        ) 

readr::write_csv(file = "fastq_paths_all.csv", x = fastq_paths.df) # write csv of all read paths with metadata

fastq_paths.df %>% # data frame of sample paths
    dplyr::filter(!grepl("Undetermined",sample_id)) %>% 
    readr::write_csv(file = "fastq_paths_samples.csv", x = .)

fastq_paths.df %>% # data frame of only Undetermined paths
    dplyr::filter(grepl("Undetermined",sample_id)) %>% 
    readr::write_csv(file = "fastq_paths_ud.csv", x = .)

# join paths df to samplesheet (drops Undetermined reads)
samdf <- dplyr::left_join(samdf, fastq_paths.df, by = c("sample_id", "fcid"))

# Add primers and target gene to sample sheet
if (stringr::str_detect(params.data_folder, "single$")) { # this is a temp fix for two datasets
    samdf <- samdf %>%
        dplyr::mutate(
            target_gene = "COI",
            pcr_primers = "fwhF2-fwhR2n",
            for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC",
            rev_primer_seq = "GTRATWGCHCCDGCTARWACWGG"
            )
} else {
    if (stringr::str_detect(params.data_folder, "dual$|full_teph$")) {
        samdf <- samdf %>%
            dplyr::mutate(
                target_gene = "COI;EIF3L",
                pcr_primers = "fwhF2-fwhR2nDac;EIF3LminiF4-EIF3lminiR4",
                for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC;GATGCGYCGTTATGCYGATGC",
                rev_primer_seq = "GTRATWGCHCCIGCTAADACHGG;TTRAAYACTTCYARATCRCC"
                )
    } else { stop("Data directory not a valid option.") }
}
### TODO: Do this step programatically, without manual input of strings

#### TODO: Import parameters from external .csv or Excel spreadsheet
#### then mutate so make sure paths are updated to absolute


# Params to add in step_add_parameters
if (stringr::str_detect(params.data_folder, "single$")) {
    params <- tibble(
        # Primer parameters
        pcr_primers = "fwhF2-fwhR2n",
        target_gene="COI",
        max_primer_mismatch=0,

        # Read filtering
        read_min_length = 20,
        read_max_length = Inf,
        read_max_ee = 1,
        read_trunc_length = 150,
        read_trim_left = 0, 
        read_trim_right = 0,
        
        # ASV filtering
        asv_min_length = 195, 
        asv_max_length = 215,
        high_sensitivity = TRUE,
        concat_unmerged = FALSE,
        genetic_code = "SGC4",
        coding = TRUE,
        phmm = paste0(projectDir,"/reference/folmer_fullength_model.rds"),
        
        # Taxonomic assignment
        idtaxa_db = paste0(projectDir,"/reference/idtaxa_bftrimmed.rds"),
        ref_fasta = paste0(projectDir,"/reference/insecta_hierarchial_bftrimmed.fa.gz"),
        idtaxa_confidence = 60,
        run_blast=TRUE,
        blast_min_identity = 97,
        blast_min_coverage = 90,
        target_kingdom = "Metazoa",
        target_phylum = "Arthropoda",
        target_class = NA,
        target_order = NA,
        target_family = NA,
        target_genus = NA,
        target_species= NA,
        
        # Sample & Taxon filtering
        min_sample_reads = 1000,
        min_taxa_reads= NA,
        min_taxa_ra = 1e-4, #1e-4 is 0.01%
            
        # General pipeline parameters
        threads = 1
    )
} else { if (stringr::str_detect(params.data_folder, "dual$|full_teph$")) {
    params <- tibble::tibble(
        # Primer parameters
        pcr_primers = c("fwhF2-fwhR2nDac", "EIF3LminiF4-EIF3lminiR4"),
        target_gene=c("COI", "EIF3L"),
        max_primer_mismatch=1,

        # Read filtering
        read_min_length = 20,
        read_max_length = Inf,
        read_max_ee = 1,
        read_trunc_length = 150,
        read_trim_left = 0,
        read_trim_right = 0,
        
        # ASV filtering
        asv_min_length = c(195, 207),
        asv_max_length = c(215, 227),
        high_sensitivity = TRUE,
        concat_unmerged = FALSE,
        genetic_code = c("SGC4", "SGC0"),
        coding = c(TRUE, TRUE),
        phmm = c(paste0(projectDir,"/reference/Bactrocera_COI.rds"), paste0(projectDir,"/reference/Bactrocera_EIF3L.rds")),
        
        # Taxonomic assignment
        idtaxa_db = c(paste0(projectDir,"/reference/COI_idtaxa.rds"),paste0(projectDir,"/reference/EIF3L_idtaxa.rds")),
        ref_fasta = c(paste0(projectDir,"/reference/COI_hierarchial.fa.gz"), paste0(projectDir,"/reference/EIF3L_hierarchial.fa.gz")),
        idtaxa_confidence = 60,
        run_blast=TRUE,
        blast_min_identity = 97,
        blast_min_coverage = 90,
        target_kingdom = "Metazoa",
        target_phylum = "Arthropoda",
        target_class = "Insecta",
        target_order = "Diptera",
        target_family = NA, 
        target_genus = NA,  
        target_species = NA,  
        
        # Sample & Taxon filtering
        min_sample_reads = c(1000, 1000),
        min_taxa_reads= NA, 
        min_taxa_ra = c(1e-4, 1e-4),
        
        # General pipeline parameters
        threads = 1
        )
    } else { stop("Data directory not a valid option.") }
}



## samdf
default_samdf <- tibble::tibble(
        sample_id = NA_character_,
        sample_name = NA_character_,
        extraction_rep = NA_integer_,
        amp_rep = NA_integer_,
        client_name = NA_character_,
        experiment_name = NA_character_,
        sample_type = NA_character_,
        collection_method = NA_character_,
        collection_location = NA_character_,
        lat_lon = NA_character_,
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
# Read in input samdf
#input_samdf <- readr::read_csv(samdf_file, show_col_types = FALSE, col_types = cols(.default = "c"))

# Make sure all columns are present and same type using a special bind operation
samdf_checked <- new_bind(default_samdf %>% filter(FALSE), samdf) 

# Check essential parameters are present
assertthat::assert_that(all(is.character(samdf_checked$sample_id)) & all(!is.na(samdf_checked$sample_id)),
                        msg = "All samples must have a sample_id in the sample_info.csv file")
assertthat::assert_that(all(is.character(samdf_checked$pcr_primers)) & all(!is.na(samdf_checked$pcr_primers)),
                        msg = "All samples must have pcr_primers in the sample_info.csv file")
assertthat::assert_that(all(is.character(samdf_checked$for_primer_seq)) & all(!is.na(samdf_checked$for_primer_seq)),
                        msg = "All samples must have a for_primer_seq in the sample_info.csv file")
assertthat::assert_that(all(is.character(samdf_checked$rev_primer_seq)) & all(!is.na(samdf_checked$rev_primer_seq)),
                        msg = "All samples must have a rev_primer_seq in the sample_info.csv file")

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
        high_sensitivity = TRUE,
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
        min_taxa_ra = 0,
        threads = 1
)
# Read in params file
#input_params <- readr::read_csv(params_file, show_col_types = FALSE, col_types = cols(.default = "c"))

# Make sure all columns are present and same type using a special bind operation
params_df <- new_bind(default_params %>% filter(FALSE), params) 

readr::write_csv(file = "params.csv", x = params_df)

# Check columns arent NA
for(i in 1:ncol(default_params)){
    param_to_check <- colnames(default_params)[i]
    if(all(is.na(params_df %>% dplyr::pull(!!param_to_check))) & !param_to_check %in% colnames(params)){
    warning(paste0("Parameter: ", param_to_check, " is NA, using default: ", default_params %>% dplyr::pull(!!param_to_check)))
    params_df <- params_df %>%
        dplyr::mutate(!!param_to_check := default_params %>% dplyr::pull(!!param_to_check))
    }
}

# Check class of all columns
for(i in 1:ncol(default_params)){
    param_to_check <- colnames(default_params)[i]
    if(!class(default_params %>% dplyr::pull(!!param_to_check)) == class(params_df %>% dplyr::pull(!!param_to_check))){
    stop(paste0("The column ", param_to_check, " in loci_params.csv file must be of class ", class(default_params %>% dplyr::pull(!!param_to_check))))
    }
}


### these checks work in a single R session but break in Nextflow

# # Check idtaxa db exists - Needs to handle multiple dbs
# check_paths <- params_df$idtaxa_db[!is.na(params_df$idtaxa_db)]%>%
#     stringr::str_split(pattern=";", n=Inf) %>% 
#     unlist()
# for(i in seq_along(check_paths)){
# assertthat::is.readable(check_paths[i])
# }

# # Check fasta exists
# check_paths <- params_df$ref_fasta[!is.na(params_df$ref_fasta)] %>%
#     stringr::str_split(pattern=";", n=Inf) %>% 
#     unlist()
# for(i in seq_along(check_paths)){
#     assertthat::is.readable(check_paths[i])
# }

# # Check phmm exists
# check_paths <- params_df$phmm[!is.na(params_df$phmm)] %>%
#     stringr::str_split(pattern=";", n=Inf) %>% 
#     unlist()
# for(i in seq_along(check_paths)){
#     assertthat::is.readable(check_paths[i])
# }


# need a way to check sequencing read inputs

## save .rds outputs
# saveRDS(object = samdf, file = "samdf.rds")
# saveRDS(object = params_df, file = "params.rds")

## split samplesheet by primer and join to params by primer
samdf_params <- samdf %>%  # split samdf loci-relevant columns across new rows, then join samdf and params  
    tidyr::separate_longer_delim(c(pcr_primers, for_primer_seq, rev_primer_seq, target_gene), delim = ";") %>% 
    dplyr::left_join(., params, by = c("pcr_primers", "target_gene"))

readr::write_csv(samdf, "samdf_original.csv") # save csv

readr::write_csv(samdf_params, "samdf_params.csv") # save csv

split_samdf <- split(samdf_params, samdf_params$pcr_primers) # split dfs by pcr_primers

for ( I in 1:length(split_samdf)) { # assign new dfs to new variables
    new_df_name <- paste0(unique(split_samdf[[I]]$pcr_primers),"_samdf")
    assign(
        paste0(unique(split_samdf[[I]]$pcr_primers),"_samdf"),
        split_samdf[[I]]
        )
    readr::write_csv( # print dfs inside work dir; maybe publish?
        x = get(new_df_name), 
        file = sprintf("%s.csv",new_df_name)
        )
}

### check parameters for validity

## check all parameters for the absence of spaces, which will affect file names later on


## BLAST database ranks
## TODO: ignore this check if 'tax_ranks' is set explicitly in loci parameters
check_param <- params_df$ref_fasta[!is.na(params_df$ref_fasta)] %>% 
    stringr::str_split(pattern=";", n=Inf) %>% 
    unlist()
for (i in seq_along(check_param)) {
    ref_fasta_path <- normalizePath(check_param[i])
    if ( stringr::str_detect(ref_fasta_path, "\\.fa\\.gz$|\\.fasta\\.gz$") ) { # if compressed
        n_ranks <- read_lines(gzfile(ref_fasta_path), n_max = 1) %>% # pull first line of .fa.gz
                stringr::str_extract(pattern = ";.*?$") %>% # extract the taxonomic rank information (including first ';')
                stringr::str_count(pattern = "[^;]+") # count the number of ranks between the ';'
    } else if ( stringr::str_detect(ref_fasta_path, "\\.fa$|\\.fasta$") ) { # if uncompressed
        n_ranks <- read_lines(ref_fasta_path, n_max = 1) %>% # pull first line of .fa.gz
                stringr::str_extract(pattern = ";.*?$") %>% # extract the taxonomic rank information (including first ';')
                stringr::str_count(pattern = "[^;]+") # count the number of ranks between the ';'
    } else { # if extension is not expected
        stop ("*** 'ref_fasta' (BLAST database) file must have '.fa(sta)' or '.fa(sta).gz' extension! ***")
    }
    # set ranks based on number of ranks (only works for 7 or 8 ranks)
    if ( n_ranks > 8 | n_ranks < 7 ) {
        stop ("*** BLAST database contains <7 or >8 ranks--please set ranks explicitly using 'params.tax_ranks'! ***")
    } else {
        message ("*** 'ref_fasta' parameter is valid ***")
    }
}

# stop(" *** stopped manually *** ") ##########################################
