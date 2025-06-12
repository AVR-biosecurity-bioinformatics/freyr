process ERROR_MODEL {
    def module_name = "error_model"
    tag "$primers; $read_group"
    label "high"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(primers), val(read_group), path(reads)


    output:   
    tuple val(direction), val(primers), val(read_group), path("*_errormodel{F,R,S}.rds"),   emit: errormodel
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
    fcid =              "${read_group}"
    pcr_primers =       "${primers}"
    reads =             "${reads}"
    threads =           "${task.cpus}"
    
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