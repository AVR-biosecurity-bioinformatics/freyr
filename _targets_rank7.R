library(targets)
library(tarchetypes)
source("R/functions.R")
options(tidyverse.quiet = TRUE)

# Load packages -----------------------------------------------------------
tar_option_set(packages = c(
  "ggplot2",
  "scales",
  "gridExtra",
  "tibble",
  "dplyr",
  "stringr",
  "tidyr",
  "purrr",
  "rlang",
  "readr",
  "magrittr",
  "patchwork",
  "vegan",
  "patchwork",
  "phyloseq",
  "DECIPHER",
  "Biostrings",
  "ShortRead",
  "savR",
  "dada2",
  "ngsReports",
  "taxreturn",
  "seqateurs"
),
imports =c(
  "taxreturn",
  "seqateurs"
),
workspace_on_error = TRUE)

# Source script files in R directory
tar_source()

# Targets pipeline
list(

  # Parameter setup ---------------------------------------------------------
  # Track input files
  tar_file(samdf_file, "sample_data/Sample_info.csv"),
  tar_file(params_file, "sample_data/loci_params.csv"),
  
  # Load input tracking files
  tar_target(samdf, {
      # define default fields
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
      input_samdf <- readr::read_csv(samdf_file, show_col_types = FALSE, col_types = cols(.default = "c"))
      
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
      return(samdf_checked)
    }),
  
  tar_target(params, {
    # define default params
    default_params <- tibble::tibble(pcr_primers = NA_character_,
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
    input_params <- readr::read_csv(params_file, show_col_types = FALSE, col_types = cols(.default = "c"))
    
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
    return(params_df)
    }, tidy_eval = FALSE
    ),
  
  # Create temporary params_primer file for tracking
  tar_file(params_primer_path,{
    out <- "output/temp/params_primer.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, max_primer_mismatch) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_primer,read_csv(params_primer_path, show_col_types = FALSE)),
  
  # Create temporary params_readfilter file for tracking
  tar_file(params_readfilter_path,{
    out <- "output/temp/params_readfilter.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, read_min_length, read_max_length, read_max_ee, 
                    read_trunc_length, read_trim_left, read_trim_right) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_readfilter, read_csv(params_readfilter_path, show_col_types = FALSE)),
  
  # Create temporary params_dada file for tracking
  tar_file(params_dada_path,{
    out <- "output/temp/params_dada.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, concat_unmerged, high_sensitivity) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_dada, read_csv(params_dada_path, show_col_types = FALSE)),
  
  # Create temporary params_asvfilter file for tracking
  tar_file(params_asvfilter_path,{
    out <- "output/temp/params_asvfilter.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, asv_min_length, asv_max_length,
                    phmm, coding, genetic_code) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_asvfilter, read_csv(params_asvfilter_path, show_col_types = FALSE)),
  
  # Create temporary params_database file for tracking
  tar_file(params_database_path,{
    out <- "output/temp/params_database.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, idtaxa_db, idtaxa_confidence, 
                    ref_fasta, blast_min_identity, blast_min_coverage, run_blast) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_database, read_csv(params_database_path, show_col_types = FALSE)),
  
  # Create temporary params_ps file for tracking
  tar_file(params_ps_path,{
    out <- "output/temp/params_ps.csv"
    params %>% 
      dplyr::select(pcr_primers, target_gene, target_kingdom, target_phylum, target_class,
                     target_order, target_family, target_genus, target_species, 
                     min_sample_reads, min_taxa_reads, min_taxa_ra) %>%
      write_csv(out)
    return(out)
  }),
  tar_target(params_ps, read_csv(params_ps_path, show_col_types = FALSE)),
  
  
  # Create directories
  tar_target(create_dirs,  step_validate_folders(getwd())),
  
  # Look for sequencing reads
  tar_files(
    fastq_path,
    purrr::map(list.dirs("data", recursive=FALSE),
                          list.files, pattern="_R[12]_", full.names = TRUE) %>%
      unlist() 
  ),

  tar_target(temp_samdf1, step_check_files(samdf, fastq_path, col_name="starting")),
  
  tar_group_by(temp_samdf1_grouped, temp_samdf1, fcid),

# Sequencing QC -------------------------------------------------------
  tar_target(seq_qc, {
             process <- temp_samdf1_grouped %>%
                dplyr::group_by(fcid) %>%
                tidyr::nest() %>%
                dplyr::mutate(seq_qc = purrr::map(fcid, step_seq_qc, quiet=FALSE, write_all=FALSE))
             out <- paste0("output/logs/", unique(process$fcid),"/",unique(process$fcid),"_flowcell_qc.pdf")
             if(is.na(process$seq_qc[[1]]$reads_pf)){
               pdf(file=out, paper="A4")
               plot.new()
               text(x=.5, y=.5, "ERROR: InterOp folder or RunInfo.xml not present") 
               try(dev.off(), silent=TRUE)
             }
             return(out)
             }, pattern = map(temp_samdf1_grouped), format="file",  iteration = "vector"),
  
  tar_target(switching_qc,{
             process <- temp_samdf1_grouped %>%
               dplyr::group_by(fcid) %>%
               tidyr::nest() %>%
               dplyr::mutate(switching_qc = purrr::map(fcid, step_switching_calc, quiet=TRUE))
             out <- paste0("output/logs/",unique(process$fcid),"/",unique(process$fcid),"_index_switching.pdf")
             if(is.na(process$switching_qc[[1]]$switch_rate)){
               pdf(file=out, paper="A4")
               plot.new()
               text(x=.5, y=.5, "ERROR: Undetermined reads file not found") 
              try(dev.off(), silent=TRUE)
             }
             return(out)
             },pattern = map(temp_samdf1_grouped), format="file", iteration = "vector"),
             

# Demultiplex and trim primers --------------------------------------------
tar_target(primer_trim,
           {
             temp_samdf1 %>%
               dplyr::left_join(params_primer, by = "pcr_primers") %>%
               dplyr::mutate(primer_trim = purrr::pmap(dplyr::select(., sample_id, for_primer_seq, rev_primer_seq, pcr_primers, fcid, max_primer_mismatch),
                                                       .f = ~step_primer_trim(sample_id = ..1, for_primer_seq=..2, rev_primer_seq=..3, pcr_primers = ..4,
                                                                              input_dir = paste0("data/",..5), output_dir =  paste0("data/",..5,"/01_trimmed"),
                                                                              qc_dir=paste0("output/logs/",..5),
                                                                              max_mismatch=..6,
                                                                              quiet = FALSE)))%>%
               dplyr::select(sample_id, sample_name, fcid, primer_trim)
           },
           pattern = map(temp_samdf1), iteration = "vector"),

# Return filepath for tracking
tar_target(primer_trim_path,
           {
             outF <- primer_trim$primer_trim %>%  
               dplyr::bind_rows()%>%
               dplyr::pull(fwd_out) 
             outR <- primer_trim$primer_trim %>%
               dplyr::bind_rows()%>%
               dplyr::pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed path
             return(c(outF,outR))
             
           },
           pattern = map(primer_trim), format="file", iteration = "vector"),

## Make temporary samdf
 tar_target(temp_samdf2, {
   temp_samdf1 %>%
                    dplyr::select(-where(is.list)) %>%
                    step_demux_samdf() %>%
                    step_check_files(primer_trim_path, col_name="trimmed")
   }),

# Filter reads ------------------------------------------------------------
  tar_target(read_filter,
             {
             temp_samdf2 %>%
                dplyr::select(sample_id, sample_name, pcr_primers, fcid) %>%
                dplyr::left_join(params_readfilter, by = "pcr_primers") %>%
                 dplyr::mutate(read_filter = purrr::pmap(dplyr::select(., sample_id, fcid, read_min_length, read_max_length, read_max_ee, 
                                                                read_trunc_length, read_trim_left, read_trim_right),
                   .f = ~step_filter_reads(
                    sample_id = ..1,
                    input_dir = paste0("data/",..2,"/01_trimmed/"),
                    output_dir = paste0("data/",..2,"/02_filtered"),
                    min_length = ..3,
                    max_length = ..4,
                    max_ee = ..5,
                    trunc_length = ..6,
                    trim_left = ..7,
                    trim_right = ..8,
                    rm.lowcomplex = 0,
                    quiet = FALSE)))%>%
                 dplyr::select(sample_id, sample_name, fcid, read_filter)
             },                  
             pattern = map(temp_samdf2), iteration = "vector"),

# Return filepath for tracking
tar_target(read_filter_path,
           {
             outF <- read_filter$read_filter %>%  
               dplyr::bind_rows()%>%
               dplyr::pull(fwd_out) 
             outR <- read_filter$read_filter %>%
               dplyr::bind_rows()%>%
               dplyr::pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed path
             return(c(outF,outR))
             
           },
           pattern = map(read_filter), format="file", iteration = "vector"),

# Make temporary samdf
 tar_target(temp_samdf3, {
   temp_samdf2 %>%
           dplyr::select(-where(is.list)) %>%
           step_check_files(read_filter_path, col_name="filtered") 
 }),




## Pre-filtering quality plots ---------------------------------------------

# Sample a random set of 5 samples for read quality plotting
#tar_target(prefilt_read_samples,{
#  group_sizes <- temp_samdf2 %>%
#    dplyr::group_by(fcid, pcr_primers) %>%
#    dplyr::group_size()
#  if(all(group_sizes > 5)){
#    n_samples <- 5
#  } else {
#    n_samples = min(group_sizes)
#  }
#  out <- temp_samdf2 %>%
#    dplyr::group_by(fcid, pcr_primers) %>%
#    dplyr::slice_sample(n=n_samples) 
#}),

tar_target(prefilt_qualplots,
           temp_samdf2 %>%
           dplyr::mutate(prefilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/01_trimmed/"),
                 truncLen=NULL, quiet = FALSE, n = 10000)
            )),
            pattern = map(temp_samdf2), iteration = "vector"),

## Write out prefilt qualplots
tar_target(write_prefilt_qualplots, {
  prefilt_qualplots %>% 
    dplyr::group_by(fcid) %>%
    tidyr::nest() %>%
    purrr::pwalk(list(fcid, data),
      .f = ~{
      pdf(file=paste0("output/logs/",..1,"/", ..1,"_prefilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
      print(..2$prefilt_qualplots)
      try(dev.off(), silent=TRUE)
      })
  out <- paste0("output/logs/",unique(prefilt_qualplots$fcid),"/", unique(prefilt_qualplots$fcid),"_prefilt_qualplots.pdf")
  return(out)
  }, format="file", iteration = "vector"),

## Post-filtering quality plots --------------------------------------------
# Get the same samples as the prefilt quality plots
#tar_target(postfilt_read_samples,{
#  out <- temp_samdf3 %>%
#    dplyr::group_by(fcid, pcr_primers) %>%
#    filter(sample_id %in% prefilt_read_samples$sample_id)
#}),

tar_target(postfilt_qualplots,
           temp_samdf3 %>%
            dplyr::mutate(postfilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/02_filtered/"), truncLen=NULL, quiet = FALSE, n = 10000)
           )),
            pattern = map(temp_samdf3), iteration = "vector"),
 
# Write out postfilt qualplots
tar_target(write_postfilt_qualplots, {
  postfilt_qualplots %>% 
    dplyr::group_by(fcid) %>%
    tidyr::nest() %>%
    purrr::pwalk(list(fcid, data),
                 .f = ~{
                   pdf(file=paste0("output/logs/",..1,"/", ..1, "_postfilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
                   print(..2$postfilt_qualplots)
                   try(dev.off(), silent=TRUE)
                 })
  out <- paste0("output/logs/",unique(postfilt_qualplots$fcid),"/", unique(postfilt_qualplots$fcid),"_postfilt_qualplots.pdf")
  return(out)
}, format="file", iteration = "vector"),

# Infer sequence variants with DADA2 --------------------------------------

  # Group temporary samdf by fcid
  tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid, pcr_primers),
  
  # Group temporary samdf by sample_ids
  tar_group_by(temp_samdf3_grouped_sample, temp_samdf3, fcid, pcr_primers, sample_id),

# Error model for forward reads
  tar_target(error_model_fwd,{
    process <- temp_samdf3_grouped %>%
      dplyr::group_by(fcid, pcr_primers) %>%
      tidyr::nest() %>%
      dplyr::mutate(error_model = purrr::pmap(dplyr::select(.,fcid, pcr_primers),
                                        .f = ~step_errormodel(fcid = ..1,
                                                         input_dir = paste0("data/",..1,"/02_filtered"),
                                                         pcr_primers = ..2,
                                                         output = paste0("output/rds/",..1,"_",..2,"_errormodelF.rds"),
                                                         qc_dir = paste0("output/logs/",..1),
                                                         read="F",
                                                         nbases=1e+08,
                                                         randomize=FALSE,
                                                         multithread=FALSE,
                                                         quiet = FALSE,
                                                         write_all = FALSE)
      ))
    return(paste0("output/rds/",unique(process$fcid),"_", unique(process$pcr_primers),"_errormodelF.rds"))
  }, format="file", pattern = map(temp_samdf3_grouped), iteration = "vector"),
 
# Error model for reverse reads
tar_target(error_model_rev,{
  process <- temp_samdf3_grouped %>%
    dplyr::group_by(fcid, pcr_primers) %>%
    tidyr::nest() %>%
    dplyr::mutate(error_model = purrr::pmap(dplyr::select(.,fcid, pcr_primers),
                                            .f = ~step_errormodel(fcid = ..1,
                                                                  input_dir = paste0("data/",..1,"/02_filtered"),
                                                                  pcr_primers = ..2,
                                                                  output = paste0("output/rds/",..1,"_",..2,"_errormodelR.rds"),
                                                                  qc_dir = paste0("output/logs/",..1),
                                                                  read="R",
                                                                  nbases=1e+08,
                                                                  randomize=FALSE,
                                                                  multithread=FALSE,
                                                                  quiet = FALSE,
                                                                  write_all = FALSE)
    ))
  return(paste0("output/rds/",unique(process$fcid),"_", unique(process$pcr_primers),"_errormodelR.rds"))
}, format="file", pattern = map(temp_samdf3_grouped), iteration = "vector"),

  # TODO:How to make it just redo one of the dadas if only one runs filtered file changed changed?

# first round denoising of forward reads
  tar_target(denoise_fwd,{
    process <- temp_samdf3_grouped_sample %>%
             dplyr::group_by(fcid, pcr_primers, sample_id) %>%
             tidyr::nest() %>%    
             dplyr::mutate(error_model = purrr::map2(fcid,pcr_primers, ~{
                readRDS(error_model_fwd[stringr::str_detect(error_model_fwd,  paste0(.x,"_",.y, "_errormodelF.rds"))]) 
              })) %>%
             dplyr::mutate(dada2 = purrr::pmap(dplyr::select(.,fcid, pcr_primers, error_model, sample_id),
                                        .f = ~step_dada2_single2(fcid = ..1,
                                                         input_dir = paste0("data/",..1,"/02_filtered"),
                                                         sample_id= ..4,
                                                         pcr_primers = ..2,
                                                         output = paste0("output/rds/", ..4,"_dada1F.rds"),
                                                         qc_dir = paste0("output/logs/",..1),
                                                         read = "F",
                                                         error_model = ..3,
                                                         multithread = FALSE,
                                                         quiet = FALSE)
             ))
        return(paste0("output/rds/", unique(process$sample_id),"_dada1F.rds"))
             }, format="file", pattern = map(temp_samdf3_grouped_sample), iteration = "vector"),

# first round denoising of reverse reads
tar_target(denoise_rev,{
  process <- temp_samdf3_grouped_sample %>%
    dplyr::group_by(fcid, pcr_primers, sample_id) %>%
    tidyr::nest() %>%    
    dplyr::mutate(error_model = purrr::map2(fcid,pcr_primers, ~{
      readRDS(error_model_rev[stringr::str_detect(error_model_rev,  paste0(.x,"_",.y, "_errormodelR.rds"))]) 
    })) %>%
    dplyr::mutate(dada2 = purrr::pmap(dplyr::select(.,fcid, pcr_primers, error_model, sample_id),
                                      .f = ~step_dada2_single2(fcid = ..1,
                                                              input_dir = paste0("data/",..1,"/02_filtered"),
                                                              sample_id= ..4,
                                                              pcr_primers = ..2,
                                                              output = paste0("output/rds/", ..4,"_dada1R.rds"),
                                                              qc_dir = paste0("output/logs/",..1),
                                                              read = "R",
                                                              error_model = ..3,
                                                              multithread = FALSE,
                                                              quiet = FALSE)
    ))
  return(paste0("output/rds/", unique(process$sample_id),"_dada1R.rds"))
}, format="file", pattern = map(temp_samdf3_grouped_sample), iteration = "vector"),

# Extract priors from forward reads
tar_target(priors_fwd,{
  process <- temp_samdf3_grouped_sample %>%
    dplyr::group_by(fcid, pcr_primers, sample_id)  %>%
    tidyr::nest() %>%    
    dplyr::mutate(priors = purrr::map(sample_id, ~{
      readRDS(denoise_fwd[stringr::str_detect(denoise_fwd,  paste0(.x,"_dada1F.rds"))])$sequence
    }))
  # Only keep the ones that appear across more than one sample
  priors <- unlist(process$priors)
  priors <- names(table(priors))[table(priors) > 1]

  # TODO: Merge with any input_priors

  saveRDS(priors, "output/rds/priorsF.rds")
  return("output/rds/priorsF.rds")
}, format="file"),

# Extract priors from reverse reads
tar_target(priors_rev,{
  process <- temp_samdf3_grouped_sample %>%
    dplyr::group_by(fcid, pcr_primers, sample_id)  %>%
    tidyr::nest() %>%    
    dplyr::mutate(priors = purrr::map(sample_id, ~{
      readRDS(denoise_rev[stringr::str_detect(denoise_rev,  paste0(.x,"_dada1R.rds"))])$sequence
    }))
  # Only keep the ones that appear across more than one sample
  priors <- unlist(process$priors)
  priors <- names(table(priors))[table(priors) > 1]
  
  # TODO: Merge with any input_priors
  
  saveRDS(priors, "output/rds/priorsR.rds")
  return("output/rds/priorsR.rds")
}, format="file"),

# Run second round of forward read denoising using priors
tar_target(denoise2_fwd,{
  process <- temp_samdf3_grouped_sample %>%
    dplyr::select(-one_of("high_sensitivity"))%>%
    dplyr::left_join(params_dada, by="pcr_primers") %>%
    dplyr::group_by(fcid, pcr_primers, sample_id, high_sensitivity) %>%
    tidyr::nest() %>%    
    dplyr::mutate(error_model = purrr::map2(fcid,pcr_primers, ~{
      readRDS(error_model_fwd[stringr::str_detect(error_model_fwd,  paste0(.x,"_",.y, "_errormodelF.rds"))]) 
    })) %>%   
    dplyr::mutate(priors = list(readRDS(priors_fwd))) %>%
    dplyr::mutate(dada2 = purrr::pmap(dplyr::select(.,fcid, pcr_primers, error_model, sample_id, priors, high_sensitivity),
                                      .f = ~{
                                        if(..6){
                                          step_dada2_single2(fcid = ..1,
                                                             input_dir = paste0("data/",..1,"/02_filtered"),
                                                             sample_id= ..4,
                                                             pcr_primers = ..2,
                                                             output = paste0("output/rds/", ..4,"_dada2F.rds"),
                                                             qc_dir = paste0("output/logs/",..1),
                                                             read = "F",
                                                             priors = ..5,
                                                             error_model = ..3,
                                                             multithread = FALSE,
                                                             quiet = FALSE) 
                                        } else {
                                          # Just copy over 
                                          saveRDS(readRDS(denoise_fwd[stringr::str_detect(denoise_fwd,  paste0("output/rds/",..4,"_dada1F.rds"))]), paste0("output/rds/",..4,"_dada2F.rds"))
                                        }
                                      }
    ))
  return(paste0("output/rds/", unique(process$sample_id),"_dada2F.rds"))
}, format="file", pattern = map(temp_samdf3_grouped_sample), iteration = "vector"),

# Run second round of reverse read denoising using priors
tar_target(denoise2_rev,{
  process <- temp_samdf3_grouped_sample %>%
    dplyr::select(-one_of("high_sensitivity"))%>%
    dplyr::left_join(params_dada, by="pcr_primers") %>%
    dplyr::group_by(fcid, pcr_primers, sample_id, high_sensitivity) %>%
    tidyr::nest() %>%    
    dplyr::mutate(error_model = purrr::map2(fcid,pcr_primers, ~{
      readRDS(error_model_rev[stringr::str_detect(error_model_rev,  paste0(.x,"_",.y, "_errormodelR.rds"))]) 
    })) %>%   
    dplyr::mutate(priors = list(readRDS(priors_rev))) %>%
    dplyr::mutate(dada2 = purrr::pmap(dplyr::select(.,fcid, pcr_primers, error_model, sample_id, priors, high_sensitivity),
                                      .f = ~{
                                        if(..6){
                                          step_dada2_single2(fcid = ..1,
                                                             input_dir = paste0("data/",..1,"/02_filtered"),
                                                             sample_id= ..4,
                                                             pcr_primers = ..2,
                                                             output = paste0("output/rds/", ..4,"_dada2R.rds"),
                                                             qc_dir = paste0("output/logs/",..1),
                                                             read = "R",
                                                             priors = ..5,
                                                             error_model = ..3,
                                                             multithread = FALSE,
                                                             quiet = FALSE) 
                                        } else {
                                          # Just copy over 
                                          saveRDS(readRDS(denoise_rev[stringr::str_detect(denoise_rev,  paste0("output/rds/",..4,"_dada1R.rds"))]), paste0("output/rds/",..4,"_dada2R.rds"))
                                        }
                                      }
    ))
  return(paste0("output/rds/", unique(process$sample_id),"_dada2R.rds"))
}, format="file", pattern = map(temp_samdf3_grouped_sample), iteration = "vector"),

# Merge reads and create seqtab by pcr_primers
tar_target(dada,{
  process <- temp_samdf3_grouped %>%
    dplyr::select(-one_of("concat_unmerged"))%>%
    dplyr::left_join(params_dada, by="pcr_primers") %>%
    dplyr::group_by(fcid, pcr_primers, concat_unmerged) %>%
    tidyr::nest() %>%
    dplyr::mutate(dada = purrr::map2(fcid, pcr_primers, ~{
      fwd_to_read <- denoise2_fwd[str_detect(denoise2_fwd, paste0(.x, ".*",.y))]
      dada_fwd <- fwd_to_read %>%
        purrr::map(function(x){readRDS(x)})
      names(dada_fwd) <- basename(fwd_to_read) %>% stringr::str_remove("_dada2F.rds")
      
      rev_to_read <- denoise2_rev[str_detect(denoise2_rev, paste0(.x, ".*",.y))]
      dada_rev <- rev_to_read %>%
        purrr::map(function(x){readRDS(x)})
      names(dada_rev) <- basename(rev_to_read) %>% stringr::str_remove("_dada2R.rds")
      return(list(dada_fwd =dada_fwd, dada_rev = dada_rev))
    })) %>%
    dplyr::mutate(dada2 = purrr::pmap(dplyr::select(.,fcid, pcr_primers, concat_unmerged, dada),
                                      .f = ~step_mergereads(fcid = ..1,
                                                       input_dir = paste0("data/",..1,"/02_filtered"),
                                                       pcr_primers = ..2,
                                                       output = paste0("output/rds/",..1,"_", ..2,"_seqtab.rds"),
                                                       qc_dir = paste0("output/logs/",..1),
                                                       quiet = FALSE,
                                                       write_all = FALSE,
                                                       concat_unmerged=..3,
                                                       dada = ..4)
    ))
},
pattern = map(temp_samdf3_grouped), iteration = "vector"),

# Return filepath for tracking
tar_target(dada_path,
           {
             return(paste0("output/rds/",unique(dada$fcid), "_", unique(dada$pcr_primers),"_seqtab.rds"))
           },
           pattern = map(dada), format="file", iteration = "vector"),

  
#  Filter ASVs -------------------------------------------------
# Filter by primer
 tar_target(filtered_seqtab, {
   temp_samdf3_grouped %>%
           dplyr::select(-one_of("asv_min_length", "asv_max_length", "phmm", "coding", "genetic_code"))%>%
           dplyr::left_join(params_asvfilter, by="pcr_primers") %>%
           dplyr::group_by(fcid, pcr_primers, asv_min_length, asv_max_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq) %>%
           tidyr::nest() %>%
           dplyr::mutate(subset_seqtab = purrr::map2(fcid, pcr_primers, ~{
              readRDS(dada_path[stringr::str_detect(dada_path,  paste0(.x, ".*",.y, "_seqtab.rds"))]) 
           })) %>%
           dplyr::ungroup()%>%
           dplyr::mutate(filtered_seqtab = purrr::pmap(dplyr::select(., fcid, pcr_primers, subset_seqtab, asv_min_length, asv_max_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq),
               .f = ~step_filter_asvs(
               seqtab = ..3,
               pcr_primers = ..2,
               output = paste0("output/rds/",..1,"_", ..2,"_seqtab.cleaned.rds"),
               qc_dir = "output/logs/",
               min_length = ..4,
               max_length = ..5,
               phmm = ..6,
               check_frame = ..7,
               genetic_code = ..8,
               primers = c(..9, ..10),
               multithread = FALSE, 
               quiet = FALSE)
         )) %>% 
         unnest_wider(filtered_seqtab) %>%
         dplyr::mutate(filtered_asvs = purrr::map(filtered_asvs, ~{
           .x %>%
             dplyr::select(-sample_id) # remove sample_id from the nested column to avoid duplicate name columns
         })) %>%
         tidyr::unnest(c(data, filtered_asvs))%>%
         dplyr::select(sample_id, sample_name, fcid, reads_starting, reads_chimerafilt, pcr_primers, reads_lengthfilt,
                       reads_phmmfilt, reads_framefilt, reads_final, plot, cleanup_summary)%>%
         dplyr::mutate(path = paste0("output/rds/",fcid, "_", pcr_primers,"_seqtab.cleaned.rds"))
 }, pattern = map(temp_samdf3_grouped), iteration = "vector"),
 
# Return filepath for tracking
tar_target(filtered_seqtab_path,
           {
             return(unique(filtered_seqtab$path))
           }, format="file"),

# Write out seqtab filtering summary csv
tar_target(write_seqtab_summary, {
  bind_rows(unique(filtered_seqtab$cleanup_summary)) %>%
    write_csv("output/logs/ASV_cleanup_summary.csv")
  out <- "output/logs/ASV_cleanup_summary.csv"
  return(out)
}, format="file", iteration = "vector"),

# Write out seqtab filtering plots
tar_target(write_seqtab_qualplots, {
  pdf("output/logs/ASV_cleanup_summary.pdf", width = 11, height = 8 , paper="a4r")
    print(unique(filtered_seqtab$plot))
  try(dev.off(), silent=TRUE)
  out <- "output/logs/ASV_cleanup_summary.pdf"
  return(out)
}, format="file", iteration = "vector"),


## Merge all loci into final seqtab ----------------------------------------
  tar_target(merged_seqtab_path, {
    process <- temp_samdf3 %>%
      tidyr::nest(data=everything()) %>%
      dplyr::mutate(final_seqtab = purrr::map(data, ~{
        seqtabs <- filtered_seqtab_path
        seqtabs <- seqtabs[seqtabs %>%
                             purrr::map_lgl(function(y){
                               any(stringr::str_detect(y, paste0(unique(.x$pcr_primers), "_seqtab.cleaned.rds")))
                             })]
        # Remove empty seqtabs
        empty <- seqtabs[purrr::map_lgl(seqtabs, function(z){
          !ncol(readRDS(z)) > 0
        })]
        if(length(empty) > 0 ){
          warning(paste0("No sequences in: ", empty, " skipping"))
        }
        if(length(seqtabs[!seqtabs %in% empty]) > 1){
          st.all <- dada2::mergeSequenceTables(tables=seqtabs[!seqtabs %in% empty])
        } else if(length(seqtabs[!seqtabs %in% empty]) == 1) {
          st.all <- readRDS(seqtabs[!seqtabs %in% empty])
        }
        saveRDS(st.all, "output/rds/seqtab_final.rds")
        return(TRUE)
      })) %>%
      tidyr::unnest(data)
    out <- paste0("output/rds/seqtab_final.rds")
    return(out)
   }, format="file", iteration = "vector"),

# Assign taxonomy ---------------------------------------------------------
  tar_file(idtaxa_db_tracked,{
           out <- params_database  %>%
             dplyr::pull(idtaxa_db) %>%
             unique()%>%
             stringr::str_split(pattern=";", n=Inf) %>% 
             unlist()
           return(out)
  }
  ),
  tar_file(ref_fasta_tracked,{
           out <- params_database  %>%
             dplyr::pull(ref_fasta) %>%
             unique() %>%
             stringr::str_split(pattern=";", n=Inf) %>% 
             unlist()
           return(out)
  }
  ),

## IDTAXA -------------------------------------------------------------------
 tar_target(tax_idtaxa,{ 
   temp_samdf3_grouped %>%
     dplyr::select(-one_of("target_gene", "idtaxa_db"))%>%
     dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene, idtaxa_db, idtaxa_confidence)) %>%
     tidyr::separate_rows(idtaxa_db, sep=";") %>%
     dplyr::group_by(fcid, target_gene, pcr_primers, idtaxa_db, idtaxa_confidence) %>%
     tidyr::nest()  %>% 
     dplyr::mutate(idtaxa_db2 = purrr::map(idtaxa_db, ~{
       idtaxa_db_tracked[stringr::str_detect(idtaxa_db_tracked, .x)]
     }))  %>%
     unnest(idtaxa_db2)%>%
     dplyr::mutate(filtered_seqtab = purrr::map2(fcid, pcr_primers, ~{
            readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path,  paste0(.x,"_", .y, "_seqtab.cleaned.rds"))]) 
     }))  %>%
     dplyr::mutate(idtaxa = purrr::pmap(list(target_gene, pcr_primers, filtered_seqtab, idtaxa_db2, idtaxa_confidence),
                                 .f = ~step_idtaxa(
                                   seqtab = ..3,
                                   database = ..4,
                                   ranks = c("Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                                   qc_dir = "output/logs/",
                                   threshold = ..5,
                                   multithread = FALSE, 
                                   quiet = FALSE,
                                   return_ids = TRUE)
     )) %>%
     mutate(idtaxa_ids = purrr::map(idtaxa,~{
       .x$ids
     })) %>%
     mutate(idtaxa = purrr::map(idtaxa,~{
       .x$tax
     })) 
   }, pattern = map(temp_samdf3_grouped), iteration = "vector"),

# Write out idtaxa objects
tar_target(idtaxa_path, {
  process <- tax_idtaxa %>%
    mutate(output = purrr::pmap(list(fcid, pcr_primers, idtaxa_db, idtaxa), ~{
      # Write out RDS of the tax table
      saveRDS(..4, paste0("output/rds/",..1,"_",..2,"_",basename(..3)  %>% stringr::str_remove("\\..*$"),"_taxtab.rds"))
    }))
  out <- unique(paste0("output/rds/",process$fcid,"_",process$pcr_primers,"_",basename(process$idtaxa_db  %>% stringr::str_remove("\\..*$")),"_taxtab.rds"))
  return(out)
}, format="file", iteration = "vector"),

## BLAST -------------------------------------------------------------------
 tar_target(tax_blast_path,
            {
           process <- temp_samdf3_grouped %>%
             dplyr::select(-one_of("target_gene", "ref_fasta"))%>%
             dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene, ref_fasta, blast_min_identity, blast_min_coverage, run_blast)) %>%
             tidyr::separate_rows(ref_fasta, sep=";") %>%
             dplyr::group_by(fcid, target_gene, pcr_primers, blast_min_identity, blast_min_coverage, ref_fasta, run_blast) %>%
             tidyr::nest()  %>% 
             dplyr::mutate(ref_fasta2 = purrr::map(ref_fasta, ~{
               ref_fasta_tracked[stringr::str_detect(ref_fasta_tracked, .x)]
             }))  %>%
             unnest(ref_fasta2) %>%
             dplyr::mutate(filtered_seqtab = purrr::map2(fcid, pcr_primers, ~{
               readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path,  paste0(.x,"_", .y, "_seqtab.cleaned.rds"))]) 
             }))  %>%
             dplyr::mutate(blast = purrr::pmap(list(fcid, target_gene, pcr_primers, filtered_seqtab, ref_fasta2, blast_min_identity, blast_min_coverage, run_blast),
                                        .f = ~{
                                        if(isTRUE(..8)){
                                        blast_res <- step_blast_tophit(
                                          seqtab = ..4,
                                          database = ..5,
                                          ranks = c("Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                                          output = paste0("output/rds/",..1,"_",..3,"_",basename(..5) %>% stringr::str_remove("\\..*$"),"_blast.rds"),
                                          qc_dir = "output/logs/",
                                          identity = ..6,  
                                          coverage=..7,
                                          evalue=1e06,
                                          max_target_seqs=5,
                                          max_hsp=5, 
                                          multithread = FALSE, 
                                          quiet = FALSE)
                                        return(blast_res)
                                        } else {
                                          blast_res <- tibble::enframe(getSequences(..4), name=NULL, value="OTU") %>%
                                            dplyr::mutate(Genus = NA_character_, Species = NA_character_) %>%
                                            column_to_rownames("OTU") %>%
                                            as.matrix()
                                          saveRDS(blast_res, paste0("output/rds/",..1,"_",..3,"_",basename(..5) %>% stringr::str_remove("\\..*$"),"_blast.rds"))
                                          return(blast_res)
                                        }
                                        }
           ))
          out <- unique(paste0("output/rds/",process$fcid, "_",process$pcr_primers, "_",basename(process$ref_fasta)  %>% stringr::str_remove("\\..*$"),"_blast.rds"))
          return(out)
        }, pattern = map(temp_samdf3_grouped), format="file", iteration = "vector"),
 
## Aggregate taxonomic assignment methods-----------------------------------------------
 tar_target(joint_tax, {
              process <- temp_samdf3_grouped %>%
                dplyr::select(-one_of("target_gene"))%>%
                dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta)) %>%
                dplyr::group_by(fcid, pcr_primers) %>%
                tidyr::nest() %>%
                dplyr::mutate(filtered_seqtab = purrr::map2(fcid, pcr_primers, ~{
                  readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path,  paste0(.x,"_", .y, "_seqtab.cleaned.rds"))]) 
                }))%>% 
                dplyr::mutate(idtaxa = purrr::pmap(list(data, fcid, pcr_primers, filtered_seqtab),
                                    .f = ~{
                                      idtaxa_dbs <- ..1 %>%
                                        tidyr::separate_rows(idtaxa_db, sep=";") %>%
                                        dplyr::pull(idtaxa_db) %>%
                                        unique() %>%
                                        basename() %>% 
                                        stringr::str_remove("\\..*$")
                                      taxtabs <- idtaxa_path[stringr::str_detect(idtaxa_path, idtaxa_dbs)& stringr::str_detect(idtaxa_path, paste0(..2, "_",..3,"_"))] %>%
                                        purrr::map(readRDS)
                                      if(length(taxtabs) == 1){
                                        out <- taxtabs[[1]]
                                      } else if(length(taxtabs) == 2){
                                        out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                                      } else if(length(taxtabs) == 3){
                                        temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                                        out <- coalesce_tax(temptax, taxtabs[[3]])
                                      }
                                      # Check that output dimensions match input
                                      if(!all(rownames(out) %in% colnames(..4))){
                                        stop("Number of ASVs classified does not match the number of input ASVs")
                                      }
                                      return(out)
                                    }),
                blast =  purrr::pmap(list(data, fcid, pcr_primers, filtered_seqtab),
                                    .f = ~{
                                      ref_fastas <- ..1 %>%
                                        tidyr::separate_rows(ref_fasta, sep=";") %>%
                                        dplyr::pull(ref_fasta) %>%
                                        unique() %>%
                                        basename() %>% 
                                        stringr::str_remove("\\..*$")
                                      taxtabs <- tax_blast_path[stringr::str_detect(tax_blast_path, ref_fastas) & stringr::str_detect(tax_blast_path, paste0(..2, "_",..3,"_"))]%>%
                                        purrr::map(readRDS)
                                      if(length(taxtabs) == 1){
                                        out <- taxtabs[[1]]
                                      } else if(length(taxtabs) == 2){
                                        out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                                      } else if(length(taxtabs) == 3){
                                        temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                                        out <- coalesce_tax(temptax, taxtabs[[3]])
                                      }
                                      # Check that output dimensions match input
                                      if(!all(rownames(out) %in% colnames(..4))){
                                        stop("Number of ASVs classified does not match the number of input ASVs")
                                      }
                                      return(out)
                                    })) %>%
                mutate(joint_tax = purrr::pmap(list(fcid, pcr_primers, idtaxa, blast, filtered_seqtab),
                                    .f = ~{
                                     out <- step_join_tax_blast(
                                       tax = ..3,
                                       blast_spp = ..4,
                                       ranks = c("Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                                       output = paste0("output/rds/",..1,"_",..2,"_taxblast.rds"),
                                       propagate_tax = TRUE)
                                     
                                     # Check that output dimensions match input
                                     if(!all(rownames(out) %in% colnames(..5))){
                                       stop("Number of ASVs classified does not match the number of input ASVs")
                                     }
                                     return(out)
                                    }
                )) %>%
                tidyr::unnest(data)
            out <- unique(paste0("output/rds/",process$fcid, "_", process$pcr_primers,"_taxblast.rds"))
            return(out)
           }, pattern = map(temp_samdf3_grouped),  format="file", iteration = "vector"),

## Merge taxonomy tables  -----------------------------------------------------
 tar_target(merged_tax,
            {
              process <- temp_samdf3 %>%
                dplyr::select(-one_of("target_gene"))%>%
                dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene)) %>%
                tidyr::nest(data=everything()) %>%
                dplyr::mutate(merged_seqtab = list(readRDS(merged_seqtab_path)))%>%
                dplyr::mutate(merged_tax = purrr::map2(data, merged_seqtab,
                    .f = ~{
                      taxtabs <- joint_tax
                      taxtabs <- taxtabs[taxtabs %>%
                                           purrr::map_lgl(function(y){
                                             any(stringr::str_detect(y, unique(paste0(.x$fcid, "_", .x$pcr_primers,"_taxblast.rds"))))
                                           })] %>%
                        purrr::map(readRDS) 
                      taxtabs <- taxtabs[sapply(taxtabs, nrow) >0 ] 

                      tax_merged <- taxtabs %>%
                        purrr::map(~{
                          .x %>%
                            tibble::as_tibble(rownames = "OTU")
                        }) %>%
                        dplyr::bind_rows() %>%
                        dplyr::distinct() # Remove any exact duplicates from save ASV being in different seqtab
                      
                      # Check for duplicated ASVs across taxtabs
                      if(any(duplicated(tax_merged$OTU))){
                        warning("Duplicated ASVs detected, selecting first occurance")
                        out <- tax_merged %>%
                          dplyr::group_by(OTU) %>%
                          dplyr::slice(1)%>%
                          tibble::column_to_rownames("OTU") %>%
                          as.matrix()
                      } else{
                        out <- tax_merged %>%
                          tibble::column_to_rownames("OTU")%>%
                          as.matrix()
                      }
                      # Check that output dimensions match input
                      if(!all(rownames(out) %in% colnames(.y))){
                        stop("Number of ASVs classified does not match the number of input ASVs")
                      }
                      saveRDS(out, "output/rds/merged_tax.rds")
                      return(out)
                    })) %>%
                tidyr::unnest(data)            
              out <- "output/rds/merged_tax.rds"
              return(out)
          }, format="file", iteration = "vector"),


## Assignment summary ---------------------------------------------------------

tar_target(assignment_plot, {
  temp_samdf3_grouped %>%
    dplyr::select(-one_of("target_gene", "idtaxa_db"))%>%
    dplyr::left_join(params_database %>% dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta)) %>%
    tidyr::separate_rows(ref_fasta, sep=";") %>%
    dplyr::group_by(fcid, target_gene, pcr_primers, idtaxa_db, ref_fasta) %>%
    tidyr::nest()  %>% 
    dplyr::mutate(ref_fasta2 = purrr::map(ref_fasta, ~{
      ref_fasta_tracked[stringr::str_detect(ref_fasta_tracked, .x)]
    }))  %>%
    unnest(ref_fasta2) %>%
    dplyr::mutate(filtered_seqtab = purrr::map2(fcid, pcr_primers, ~{
      #TODO this may need to be merged by PCR primer
      readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path,  paste0(.x,"_", .y, "_seqtab.cleaned.rds"))]) 
    }))%>% 
    dplyr::mutate(tax = purrr::map2(fcid,pcr_primers, ~{
      readRDS(joint_tax[stringr::str_detect(joint_tax, paste0(.x,"_",.y,"_taxblast.rds"))])%>% 
        seqateurs::unclassified_to_na(rownames=FALSE) %>%
        dplyr::mutate(lowest = seqateurs::lowest_classified(.)) 
    })) %>%
    dplyr::mutate(blast = purrr::pmap(list(pcr_primers, filtered_seqtab, ref_fasta2),
                               .f = ~{
                                 seqmap <- tibble::enframe(getSequences(..2), name = NULL, value="OTU") %>%
                                   mutate(name = paste0("SV", seq(length(getSequences(..2)))))
                                 seqs <- taxreturn::char2DNAbin(seqmap$OTU)
                                 names(seqs) <- seqmap$name
                                 
                                 if(length(seqs) > 0){
                                   blast_top_hit(
                                   query = seqs,
                                   db = ..3,
                                   identity=60,
                                   coverage=80,
                                   ranks = c("Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")) %>% 
                                     mutate(blastspp = paste0(Genus, " ", Species)) %>%
                                     dplyr::select(name = qseqid, acc, blastspp, pident, total_score, max_score, evalue, qcovs) %>%
                                     left_join(seqmap) %>%
                                     dplyr::select(-name)
                                 } else {
                                   tibble::enframe(seqs, name=NULL, value="OTU") %>%
                                     dplyr::mutate(acc = NA_character_,
                                                   blastspp = NA_character_,
                                                   pident = NA_real_, 
                                                   length = NA_integer_,
                                                   evalue = NA_real_,
                                                   qcovs = NA_integer_) 
                                 }
    })) %>%
    dplyr::mutate(joint = purrr::pmap(list(blast, tax),
                               .f = ~{
                                 if(nrow(..1) > 0 & nrow(..2) > 0){
                                 ..1 %>%
                                   dplyr::left_join(..2, by="OTU")
                                 } else {
                                  NULL
                                 }
    })) %>%
    dplyr::mutate(plot = purrr::pmap(list(fcid, pcr_primers, joint, idtaxa_db, ref_fasta),
                              .f= ~{
                                # ADD TITLES HERE!
                                if(!is.null(..3)){
                                  
                                  cols <- c(Root = "#D53E4F",
                                            Kingdom = "#F46D43",
                                            Phylum = "#FDAE61",
                                            Class = "#FEE08B",
                                            Order = "#E6F598",
                                            Family = "#ABDDA4",
                                            Genus = "#66C2A5",
                                            Species = "#3288BD") 
                                  ..3 %>%
                                    dplyr::select(pident, rank = lowest) %>%
                                    dplyr::mutate(rank = factor(rank, levels = c("Root","Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
                                    ggplot(aes(x=pident, fill=rank))+ 
                                    geom_histogram(colour="black", binwidth = 1, position = "stack") + 
                                    labs(title = paste0(..1, "  ", ..2, " Top hit identity distribution"),
                                         subtitle = paste0("IDTAXA database:", ..4, " BLAST database:", ..5),
                                         x = "BLAST top hit % identity",
                                         y = "Sequence Variants") + 
                                    scale_x_continuous(breaks=seq(60,100,2)) +
                                    scale_fill_manual(name = "Taxonomic \nAssignment", values = cols)+
                                    theme_bw()+
                                    theme(
                                      strip.background = element_rect(colour = "black", fill = "lightgray"),
                                      strip.text = element_text(size=9, family = ""),
                                      plot.background = element_blank(),
                                      text = element_text(size=9, family = ""),
                                      axis.text = element_text(size=8, family = ""),
                                      legend.position = "right",
                                      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                                      panel.grid = element_line(size = rel(0.5)),
                                    ) 
                                  } else {
                                    NULL
                                  }
                              }))
}, pattern = map(temp_samdf3_grouped), iteration = "vector"),

# Write out assignment plot
tar_target(write_assignment_plot, {
  out <- paste0("output/logs/taxonomic_assignment_summary.pdf")
  if(!all(sapply(assignment_plot$plot, is.null))){
    pdf(file=out, width = 11, height = 8 , paper="a4r")
      print(unique(assignment_plot$plot))
    try(dev.off(), silent=TRUE)
  } else{
    pdf(file=out, width = 11, height = 8 , paper="a4r")
      plot.new()
      text(x=.5, y=.5, "ERROR: No blast hits to reference fasta - assignment plot not created") 
    try(dev.off(), silent=TRUE)
  }
  return(out)
}, format="file", iteration = "vector"),


# Write out assignment results
tar_target(tax_summary, {
  idtaxa_summary <- tax_idtaxa %>%
    ungroup()%>%
    mutate(summary = purrr::map2(idtaxa_ids, idtaxa, ~{
      .x %>%
      purrr::map_dfr(function(x){
          taxa <- paste0(x$taxon,"_", x$confidence)
          taxa[startsWith(taxa, "unclassified_")] <- NA
          data.frame(t(taxa)) %>%
            magrittr::set_colnames(c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species")[1:ncol(.)])
        }) %>%
        mutate_all(function(y){
          name <- y %>%
            str_remove("_[0-9].*$")
          conf <- y %>%
            str_remove("^.*_") %>%
            str_trunc(width=6, side="right", ellipsis = "")
          paste0(name, "_", conf, "%") 
        }) %>%
        mutate_all(~ na_if(., "NA_NA%"))  %>%
        mutate(OTU = rownames(.y))
    })) %>%
    dplyr::select(pcr_primers,idtaxa_db, summary)%>%
    unnest(summary) %>%
    distinct()
  
  if(!any(sapply(assignment_plot$joint, is.null))){
    blast_summary <- assignment_plot %>%
      dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta, joint)%>% 
      unnest(joint)  %>%
      dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta, OTU, acc, blast_top_hit = blastspp, 
                    blast_identity = pident, blast_evalue = evalue, blast_total_score = total_score, blast_max_score = max_score, blast_qcov = qcovs)
  } else {
    blast_summary <- assignment_plot %>%
      dplyr::select(pcr_primers, target_gene, idtaxa_db, ref_fasta, joint)%>% 
      unnest(joint)
  }
  
  summary_table <- idtaxa_summary %>%
    left_join(blast_summary) %>%
    dplyr::select(any_of(c("OTU", "pcr_primers", "target_gene", "idtaxa_db", "ref_fasta", 
                  "Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species",
                  "blast_top_hit", "blast_identity", "blast_qcov","blast_evalue"
    )))
  
  out <- paste0("output/logs/taxonomic_assignment_summary.csv")
  write_csv(summary_table, out)
  return(out)
}, format="file", iteration = "vector"),

# Create phyloseq object --------------------------------------------------
  tar_target(ps,{
    process <- step_phyloseq(
       seqtab = merged_seqtab_path,
       taxtab = merged_tax,
       samdf = temp_samdf2 %>% dplyr::select(-any_of(c("starting_fwd", "starting_rev", 
                                                     "trimmed_fwd",	"trimmed_rev"))),
       seqs=NULL,
       phylo=NULL,
       name_variants=TRUE)
    out <- "output/rds/ps.rds"
    saveRDS(process, out)
    return(out)
  }, format="file", iteration = "vector"),

## Output unfiltered results -----------------------------------------------
  tar_target(ps_summary, {
    ps_obj <- ps %>%
      readRDS()
    
    out <- c(step_output_summary(ps_obj, out_dir="output/results/unfiltered", type="unfiltered"),
             step_output_ps(ps_obj, out_dir="output/results/unfiltered", type="unfiltered"))
    return(out)
  }, format="file", iteration = "vector"),



# Accumulation curve ------------------------------------------------------
tar_target(accumulation_curve, {
  ps_obj <- ps %>%
    readRDS()
  gg.acc_curve <- rareplot(ps_obj, step="auto", threshold = max(params_ps$min_sample_reads))
  out <- "output/logs/accumulation_curve.pdf"
  pdf(file=out, width = 11, height = 8 , paper="a4r")
    print(gg.acc_curve)
  try(dev.off(), silent=TRUE)
  return(out)
}, format="file", iteration = "vector"),

# Filter phyloseq ---------------------------------------------------------
# Taxonomic and minimum abundance filtering
tar_target(ps_filtered,{
     # if multiple primers were used - split ps into different primers
    process <- params_ps %>%
       dplyr::mutate(ps_obj = purrr::map(pcr_primers,
              .f = ~{
                    physeq <- ps %>%
                      readRDS()
                    new_samdat <- as(phyloseq::sample_data(physeq), "data.frame") %>%
                      dplyr::filter(pcr_primers == .x)
                    if(nrow(new_samdat) > 0){
                      phyloseq::sample_data(physeq) <- phyloseq::sample_data(new_samdat)
                      physeq <- physeq %>% 
                        filter_taxa(function(x) mean(x) > 0, TRUE) #Drop missing taxa from table
                      return(physeq)
                    } else{
                      return(NULL)
                    }
                      })) %>%
    dplyr::mutate(ps_filt = purrr::pmap(dplyr::select(., ps_obj, target_kingdom, target_phylum, target_class,
                                               target_order, target_family, target_genus, target_species, 
                                               min_sample_reads, min_taxa_reads, min_taxa_ra),
            .f = ~{
              if(!is.null(..1)){
                ..1 %>%
                step_filter_phyloseq(
                  kingdom = ..2,
                  phylum = ..3,
                  class = ..4,
                  order = ..5,
                  family = ..6,
                  genus = ..7,
                  species = ..8,
                  min_sample_reads=..9, 
                  min_taxa_reads = ..10,
                  min_taxa_ra = ..11,
                  quiet=FALSE
                )
              }else {
                return(NULL)
              }
              }))
    ps_merged <- merge_phyloseq_new(process$ps_filt)
    out <- "output/rds/ps_filtered.rds"
    saveRDS(ps_merged, out)
    return(out)
  }, format="file", iteration = "vector"),

## Output filtered results -------------------------------------------------
tar_target(ps_filt_summary, {
  ps_obj <- ps_filtered %>%
    readRDS()
  
  out <- c(step_output_summary(ps_obj, out_dir="output/results/filtered", type="filtered"),
           step_output_ps(ps_obj, out_dir="output/results/filtered", type="filtered"))
  return(out)
}, format="file", iteration = "vector"),

# Read tracking ------------------------------------------------------
tar_target(read_tracking, {
  
  # Unfiltered phyloseq
  ps_obj <- readRDS(ps)
  phyloseq::tax_table(ps_obj) <- phyloseq::tax_table(ps_obj)  %>%
    as("matrix") %>% 
    as.data.frame() %>%
    seqateurs::unclassified_to_na() %>%
    as.matrix()

  # filtered phyloseq
  ps_filt_obj <- readRDS(ps_filtered)
  phyloseq::tax_table(ps_filt_obj) <- phyloseq::tax_table(ps_filt_obj)%>%
    as("matrix") %>% 
    as.data.frame() %>%
    seqateurs::unclassified_to_na() %>%
    as.matrix()
      
  read_tracker <- temp_samdf2 %>%
    dplyr::select(sample_name, sample_id, fcid, pcr_primers) %>%
    distinct()%>%
    dplyr::left_join(primer_trim %>%
                tidyr::unnest(primer_trim) %>%
                dplyr::select(sample_name, fcid, trimmed_input, trimmed_output, fwd_out) %>%
                dplyr::mutate(sample_id = basename(fwd_out) %>% stringr::str_remove("_S[0-9]+_R[1-2]_.*$")) %>%
                dplyr::select(sample_id, fcid, input_reads=trimmed_input, trimmed=trimmed_output), 
                by = c("sample_id", "fcid")) %>%
    dplyr::left_join(read_filter %>%
                tidyr::unnest(read_filter) %>%
                dplyr::select(sample_id, fcid, filtered = filter_output),
              by = c("fcid", "sample_id")) %>%
    dplyr::left_join(dada %>% 
                ungroup()%>%
                tidyr::unnest(dada2) %>% 
                dplyr::select(fcid, sample_id, denoised=merged),
              by = c("fcid", "sample_id")
              ) %>%
    dplyr::left_join(filtered_seqtab %>% 
                dplyr::select(sample_id, fcid, chimerafilt=reads_chimerafilt,
                              lengthfilt= reads_lengthfilt, phmmfilt=reads_phmmfilt, framefilt = reads_framefilt),
              by = c("fcid", "sample_id")
              ) %>%
    dplyr::left_join(psmelt(ps_obj) %>% # Could replace this with seqateurs::lowest_classified?
                dplyr::filter(Abundance > 0) %>%
                dplyr::select(sample_id, fcid, Abundance, any_of(colnames(phyloseq::tax_table(ps_obj)))) %>%
                pivot_longer(cols=any_of(colnames(phyloseq::tax_table(ps_obj))), 
                             names_to = "rank",
                             values_to="name") %>%
                filter(!is.na(name))%>%
                dplyr::group_by(sample_id, rank) %>%
                summarise(Abundance = sum(Abundance)) %>%
                pivot_wider(names_from="rank",
                            values_from="Abundance")%>%
                rename_with(~stringr::str_to_lower(.), everything()) %>%
                rename_with(~stringr::str_c("classified_", .), -sample_id), by="sample_id")  %>%
    dplyr::left_join(psmelt(ps_filt_obj) %>% # Could replace this with seqateurs::lowest_classified?
                       dplyr::filter(Abundance > 0) %>%
                       dplyr::group_by(sample_id) %>%
                       summarise(sample_taxon_filt = sum(Abundance)), by="sample_id")  %>%
    dplyr::select(any_of(c(
      "sample_name","sample_id", "pcr_primers", "fcid", "input_reads", "trimmed", "filtered",
      "denoised", "chimerafilt", "lengthfilt", "phmmfilt", "framefilt", 
      "classified_root", "classified_kingdom", "classified_phylum","classified_class",
      "classified_order", "classified_family", "classified_genus", "classified_species", "sample_taxon_filt"
    )))
  
  write_csv(read_tracker, "output/logs/read_tracker.csv")

  gg.read_tracker <- read_tracker %>%
    pivot_longer(cols = -c("sample_name","sample_id", "fcid", "pcr_primers"),
                 names_to = "step",
                 values_to="reads") %>%
    group_by(sample_name) %>%
    group_modify(~{
      # When a sample name shares multiple sample ids, select a single sample to avoid double counting
      if(length(unique(.x$pcr_primers)) >1){
        .x %>%
          dplyr::filter(!step == "input_reads") %>%
          bind_rows(.x %>%
                      dplyr::filter(step == "input_reads") %>%
                      mutate(pcr_primers="Mixed",
                             sample_id=NA_character_) %>%
                      dplyr::slice(1))
      } else {
        .x
      }
    }) %>%
    dplyr::mutate(step = factor(step, levels=c(
      "input_reads", "trimmed", "filtered",
      "denoised", "chimerafilt", "lengthfilt", "phmmfilt", "framefilt", 
      "classified_root", "classified_kingdom", "classified_phylum","classified_class",
      "classified_order", "classified_family", "classified_genus", "classified_species", "sample_taxon_filt"
    ))) %>%
    ggplot(aes(x = step, y = reads, fill=pcr_primers))+
    geom_col() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))+
    facet_grid(fcid~.)+
    theme_bw()+
    theme(
      strip.background = element_rect(colour = "black", fill = "lightgray"),
      strip.text = element_text(size=9, family = ""),
      axis.text.x =element_text(angle=45, hjust=1, vjust=1),
      plot.background = element_blank(),
      text = element_text(size=9, family = ""),
      axis.text = element_text(size=8, family = ""),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      panel.grid = element_line(linewidth = rel(0.5)),
    ) +
    labs(x = "Pipeline step",
         y = "Reads retained",
         fill = "PCR primers")
  pdf(file="output/logs/read_tracker.pdf", width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker)
  try(dev.off(), silent=TRUE)
  
  return(c("output/logs/read_tracker.csv", "output/logs/read_tracker.pdf"))
}, format="file", iteration = "vector")

)

