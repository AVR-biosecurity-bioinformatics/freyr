process PHYLOSEQ_FILTER {
    def module_name = "phyloseq_filter"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(ps), path(filters_tibble), val(process_params)

    output:
    path("*.csv"),                                                          emit: csvs
    tuple val(primers), path("ps_filtered_*.rds"),                          emit: ps 
    path("*.fasta"),                                                        emit: asv_fasta

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primers =               "${primers}"
    ps =                    "${ps}"
    filters_tibble =        "${filters_tibble}"
    target_kingdom =        "${process_params.target_kingdom}"
    target_phylum =         "${process_params.target_phylum}"
    target_class =          "${process_params.target_class}"
    target_order =          "${process_params.target_order}"
    target_family =         "${process_params.target_family}"
    target_genus =          "${process_params.target_genus}"
    target_species =        "${process_params.target_species}"
    min_sample_reads =      "${process_params.min_sample_reads}"
    min_taxa_reads =        "${process_params.min_taxa_reads}"
    min_taxa_ra =           "${process_params.min_taxa_ra}"
    
    ## global variables
    projectDir = "$projectDir"
    params_dict = "$params"
    
    tryCatch({
    ### source functions and themes, load packages, and import Nextflow params
    ### from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
    )
    }, finally = {
    ### save R environment for debugging
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """

}