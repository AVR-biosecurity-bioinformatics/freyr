process PARSE_INPUTS {
    def module_name = "parse_inputs"
    tag "Whole dataset"
    // label 

    input:
    val(samplesheet)
    val(loci_params)

    output: 

    publishDir "${projectDir}/output/modules/${module_name}"

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    samplesheet =           "${samplesheet}"
    loci_params =           "${loci_params}"

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