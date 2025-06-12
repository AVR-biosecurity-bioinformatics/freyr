process DENOISE {
    def module_name = "denoise"
    tag "$sample_primers"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path(errormodel), path(priors)
    val(n_pass)

    output:
    tuple val(direction), val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path("*_dada{1,2}{F,R,S}.rds"), emit: seq

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    def priors_path = priors.name != "NO_FILE" ? priors.name : "NO_FILE"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    direction =         "${direction}"
    fcid =              "${read_group}"
    pcr_primers =       "${primers}"
    sample_id =         "${sample_primers}"
    reads =             "${reads}"
    errormodel =        "${errormodel}"
    n_pass =            "${n_pass}"
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