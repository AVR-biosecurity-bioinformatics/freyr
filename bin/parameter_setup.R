#!/usr/bin/env Rscript

### need to make sure slashes are added correctly after data_dir
# find all directories within data folder
if (!exists("data_dir")) {data_dir="data"} # if data_dir not defined, use "data/"
print(paste0("data_dir = ", data_dir))
data_dir_abs <- paste0(projectDir,"/",data_dir)
print(paste0("data_dir_abs = ", data_dir_abs))
runs <- dir(data_dir_abs) # define data directory in module .nf file
SampleSheet <- list.files(paste0(data_dir_abs,"/", runs), pattern= "SampleSheet", full.names = TRUE) # list all sampleshets
runParameters <- list.files(paste0(data_dir_abs,"/", runs), pattern= "[Rr]unParameters.xml", full.names = TRUE) # list all run parameter files

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4") %>%
distinct()

# Check that sample_ids contain fcid, if not; attach
samdf <- samdf %>%
mutate(sample_id = case_when(
    !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
    TRUE ~ sample_id
))

# Check that samples match samplesheet
fastqFs <- 
    purrr::map(list.dirs(data_dir_abs, recursive=FALSE),
                    list.files, pattern="_R1_", full.names = TRUE) %>%
    unlist() %>%
    str_remove(pattern = "^(.*)\\/") %>%
    str_remove(pattern = "(?:.(?!_S))+$")

# Filter undetermined reads from sample sheet
fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

# write_csv(as.data.frame(fastqFs), "fastqFs.csv")

# Check for fastq files that are missing from samplesheet
if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") }

# Check for sample_ids that dont have a corresponding fastq file
if (length(setdiff(samdf$sample_id, fastqFs)) > 0) {
warning(paste0("The fastq file: ",
                setdiff(samdf$sample_id, fastqFs),
                " is missing, dropping from samplesheet \n")) 
samdf <- samdf %>%
    filter(!sample_id %in% setdiff(samdf$sample_id, fastqFs))
}

### TODO: Add .fastq paths to the sample sheet in additional fwd and rev columns
fastq_paths <- 
    list.files(data_dir_abs, pattern = "fastq.gz$|fq.gz$", recursive = T, full.names = T)

fastq_paths.fwd <- 
    fastq_paths %>% stringr::str_extract("_R1_") %>%
    as.data.frame() %>%
    write_delim("fwd.txt")

quit(status = 0)

# Add primers to sample sheet
samdf <- samdf %>%
    mutate(pcr_primers = "fwhF2-fwhR2n",
            for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC",
            rev_primer_seq = "GTRATWGCHCCDGCTARWACWGG"
            )

# write_csv(samdf, paste0(projectDir, "/sample_data/Sample_info.csv"))


# Params to add in step_add_parameters
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

# write_csv(params, "sample_data/loci_params.csv")




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

################ Creating temp files needed during targets but not now

# Create params_primer file for later
params %>% 
    dplyr::select(pcr_primers, target_gene, max_primer_mismatch) %>%
    write_csv("params_primer.csv")

# Create params_readfilter file for later
params %>% 
    dplyr::select(pcr_primers, target_gene, read_min_length, read_max_length, read_max_ee, 
                read_trunc_length, read_trim_left, read_trim_right) %>%
    write_csv("params_readfilter.csv")

# Create params_dada file for later
params %>% 
    dplyr::select(pcr_primers, target_gene, concat_unmerged, high_sensitivity) %>%
    write_csv("params_dada.csv")

# Create params_asvfilter file for later
params %>% 
    dplyr::select(pcr_primers, target_gene, asv_min_length, asv_max_length,
                phmm, coding, genetic_code) %>%
    write_csv("params_asvfilter.csv")

# Create temporary params_database file for tracking
params %>% 
    dplyr::select(pcr_primers, target_gene, idtaxa_db, idtaxa_confidence, 
                ref_fasta, blast_min_identity, blast_min_coverage, run_blast) %>%
    write_csv("params_database.csv")

# Create temporary params_ps file for tracking
params %>% 
    dplyr::select(pcr_primers, target_gene, target_kingdom, target_phylum, target_class,
                    target_order, target_family, target_genus, target_species, 
                    min_sample_reads, min_taxa_reads, min_taxa_ra) %>%
    write_csv("params_ps.csv")

# need a way to check sequencing read inputs

## save .rds outputs
saveRDS(object = samdf, file = "samdf.rds")
saveRDS(object = params, file = "params.rds")

write_csv(file = "samdf.csv", x = samdf)
write_csv(file = "params.csv", x = params)




