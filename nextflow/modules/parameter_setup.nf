process PARAMETER_SETUP {
    // Top is used for process directives like tag, label, time etc. 
    tag //$some_var
    /* Tags are used to associate each running process with a custom label, perhaps based on the file being used as input. 
    See here: https://www.nextflow.io/docs/latest/process.html#tag */
    
    // label 'process_low'
    /* Labels can be used to group processes based on certain qualities.
    This is commonly used to group by resource allocation and nf-core has standard ones. */

    input:



    output:
    /* using ", emit: name" at the end of an output line can be used to name this output later.
    For example, in ampliseq, RENAME_RAW_DATA_FILES uses ", emit: fastq" for the renamed files.
    In the main workflow, when RENAME_RAW_DATA_FILES is called, the fastq output can be used by another process
    as "RENAME_RAW_DATA_FILES.out.fastq". */

    /* The ampliseq modules/processes also output the version of the software used with additional
    bit of code. I could do this later but not worth trying for first implementation. 
    */

    path "sample_data/loci_params.csv",     emit: loci_params
    path "sample_data/Sample_info.csv",     emit: sample_info


    when: 

    /* The when block allows for conditions to be set on the running of the process, for example
    if the input file name or other characteristic matches a pattern. 
    HOWEVER: this is usually better done in the workflow block (ie. in the code where the process
    is called) so the module is more portable. But there could be good uses for it.
    The ampliseq pipeline uses "task.ext.when == null || task.ext.when" in its when blocks, but
    I'm not sure what they mean yet. 
    */


    script:
    /* This is an R script, so needs to start with "#!/usr/bin/env Rscript".
    If I want to source the "functions.R" file every time, I can keep the path to that file as 
    a variable, perhaps "r_functions", and define that in the main pipeline. */ 
    
    /* This is based on the "prepare_inputs.R" script in the R pipeline. Note that "$" in R needs to be escaped like "\$" to avoid 
    it being parsed as a Nextflow variable. */
    
    """
    #!/usr/bin/env Rscript

    #set working directory to projectDir so everything works
    setwd("${projectDir}")

    # load R libraries
    source("running_scripts/_targets_packages.R") # make sure this is the current location
    # also check that I can't just load a few packages like purrr, dplyr and tibble to save loading all the packages in the targets pipeline
    source("${functions_r}") # this should source the "functions.R" script

    # input data
    runs <- dir("data/") #Find all directories within data
    SampleSheet <- list.files(paste0("data/", runs), pattern= "SampleSheet", full.names = TRUE)
    runParameters <- list.files(paste0("data/", runs), pattern= "[Rr]unParameters.xml", full.names = TRUE)

    # Create samplesheet containing samples and run parameters for all runs
    samdf <- create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4") %>%
    distinct()

    # Check that sample_ids contain fcid, if not; attatch
    samdf <- samdf %>%
    mutate(sample_id = case_when(
        !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
        TRUE ~ sample_id
    ))

    # Check that samples match samplesheet
    fastqFs <- 
        purrr::map(list.dirs("data", recursive=FALSE),
                        list.files, pattern="_R1_", full.names = TRUE) %>%
        unlist() %>%
        str_remove(pattern = "^(.*)\\/") %>%
        str_remove(pattern = "(?:.(?!_S))+$")

    # Filter undetermined reads from sample sheet
    fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

    # Check for fastq files that are missing from samplesheet
    if (length(setdiff(fastqFs, samdf\$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf\$sample_id), " are not in the sample sheet") }

    # Check for sample_ids that dont have a corresponding fastq file
    if (length(setdiff(samdf\$sample_id, fastqFs)) > 0) {
    warning(paste0("The fastq file: ",
                    setdiff(samdf\$sample_id, fastqFs),
                    " is missing, dropping from samplesheet \n")) 
    samdf <- samdf %>%
        filter(!sample_id %in% setdiff(samdf\$sample_id, fastqFs))
    }

    # Write out sample tracking sheet
    write_csv(samdf, "sample_data/Sample_info.csv")

    # Add primers to sample sheet
    samdf <- samdf %>%
        mutate(pcr_primers = "fwhF2-fwhR2n",
            for_primer_seq = "GGDACWGGWTGAACWGTWTAYCCHCC",
            rev_primer_seq = "GTRATWGCHCCDGCTARWACWGG"
            )

    write_csv(samdf, "sample_data/Sample_info.csv")


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
        phmm = "reference/folmer_fullength_model.rds",
        
        # Taxonomic assignment
        idtaxa_db = "reference/idtaxa_bftrimmed.rds",
        ref_fasta = "reference/insecta_hierarchial_bftrimmed.fa.gz",
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

    write_csv(params, "sample_data/loci_params.csv")

    """

}