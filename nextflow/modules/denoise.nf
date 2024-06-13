process DENOISE {
    def module_name = "denoise"
    tag "$pcr_primers; $meta.sample_id"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(direction), val(pcr_primers), val(fcid), val(meta), path(reads), path(errormodel), path(priors)
    val(n_pass)

    output:
    tuple val(direction), val(pcr_primers), val(fcid), val(meta), path(reads), path("*_dada{1,2}{F,R}.rds"), emit: seq

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
    fcid =              "${fcid}"
    pcr_primers =       "${pcr_primers}"
    sample_id =         "${meta.sample_id}"
    reads =             "${reads}"
    errormodel =        "${errormodel}"
    n_pass =            "${n_pass}"
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