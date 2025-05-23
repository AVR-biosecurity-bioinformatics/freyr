process PRIORS {
    def module_name = "priors"
    tag "$pcr_primers; $fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(pcr_primers), val(fcid), path(priors)

    output:   
    tuple val(direction), val(pcr_primers), val(fcid), path("*_priors{F,R,S}.rds"),              emit: priors

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    direction =         "${direction}"
    fcid =              "${fcid}"
    pcr_primers =       "${pcr_primers}"
    priors =            "${priors}"
    
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