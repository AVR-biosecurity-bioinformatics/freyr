#!/usr/bin/env Rscript

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
samdf_checked <- new_bind(default_samdf %>% filter(FALSE), input_samdf) 

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
params_df <- new_bind(default_params %>% filter(FALSE), input_params) 

# Check columns arent NA
for(i in 1:ncol(default_params)){
    param_to_check <- colnames(default_params)[i]
    if(all(is.na(params_df %>% dplyr::pull(!!param_to_check))) & !param_to_check %in% colnames(input_params)){
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

# Check idtaxa db exists - Needs to handle multiple dbs
check_paths <- params_df$idtaxa_db[!is.na(params_df$idtaxa_db)]%>%
    stringr::str_split(pattern=";", n=Inf) %>% 
    unlist()
for(i in seq_along(check_paths)){
assertthat::is.readable(check_paths[i])
}

# Check fasta exists
check_paths <- params_df$ref_fasta[!is.na(params_df$ref_fasta)] %>%
    stringr::str_split(pattern=";", n=Inf) %>% 
    unlist()
for(i in seq_along(check_paths)){
    assertthat::is.readable(check_paths[i])
}

# Check phmm exists
check_paths <- params_df$phmm[!is.na(params_df$phmm)] %>%
    stringr::str_split(pattern=";", n=Inf) %>% 
    unlist()
for(i in seq_along(check_paths)){
    assertthat::is.readable(check_paths[i])
}


## outputs

write.csv("input_samdf.csv", input_samdf)
write.csv("params_df.csv", params_df)
