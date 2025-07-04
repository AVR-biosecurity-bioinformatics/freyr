process ACCUMULATION_CURVE {
    def module_name = "accumulation_curve"
    tag "$primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(ps_file), path(filters_tibble), val(process_params)

    output:
    path("accumulation_curve_*.pdf"),                                          emit: pdf

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primers =               "${primers}"
    ps_file =               "${ps_file}"
    min_sample_reads =      "${process_params.min_sample_reads}"
    
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