process ERROR_MODEL {
    def module_name = "error_model"
    tag "$pcr_primers; $fcid"
    // label:  

    input:
    tuple val(direction), val(pcr_primers), val(fcid), path(reads)


    output:   
    tuple val(direction), val(pcr_primers), val(fcid), path("*_errormodel{F,R}.rds"),   emit: errormodel
    path("*_errormodel.pdf"),                                           emit: plot

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
    reads =             "${reads}"
    
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