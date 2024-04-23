process DADA_PRIORS {
    def module_name = "dada_priors"
    tag "$fcid; $pcr_primers"
    // label:  

    input:
    tuple val(direction), val(fcid), val(pcr_primers), path(priors)

    output:   
    tuple val(direction), val(fcid), val(pcr_primers), path("*_priors{F,R}.rds"),              emit: priors

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
    
    ### source functions and themes, load packages, and import Nextflow params
    ### from "bin/process_start.R"
    sys.source("${projectDir}/bin/process_start.R", envir = .GlobalEnv)

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
    )
    ### save .RData for debugging
    if ("${params.rdata}" == "true") {
        save.image()
    } else {
        NULL
    }

    """

}