process MERGE_TAX {
    def module_name = "merge_tax"
    tag "Whole dataset"
    // label:  
    cache 'deep'

    input:
    tuple val(fcid), val(pcr_primers), val(meta), path(taxtab)

    output:
    tuple val(fcid), val(pcr_primers), val(meta), path("merged_tax.rds")
    path("*.rds")

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid =                  "${fcid}"
    pcr_primers =           "${pcr_primers}"
    meta =                  "${meta}"
    taxtab =                "${taxtab}"
    
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