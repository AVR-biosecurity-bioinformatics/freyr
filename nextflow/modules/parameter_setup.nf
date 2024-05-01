process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    tag "Whole dataset"
    // label: 

    input:


    output:  
    path("params.csv"),                emit: loci_params
    path("samdf_original.csv"),        emit: samdf  
    path("samdf_*.csv"),               emit: samdf_locus
    path("samdf_params.csv"),           emit: samdf_params

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables

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