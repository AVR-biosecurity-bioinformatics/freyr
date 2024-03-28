process DENOISE {
    def module_name = "denoise"
    tag "$fcid; $pcr_primers; $meta.sample_id"
    // label:  

    input:
    tuple val(direction), val(fcid), val(pcr_primers), val(meta), path(reads), path(errormodel)

    output:   


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
    sample_id =         "${meta.sample_id}"
    reads =             "${reads}"
    errormodel =        "${errormodel}"
    
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

    """

}