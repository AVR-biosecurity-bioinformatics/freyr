process SPLIT_LOCI {
    def module_name = "split_loci"
    // tag: 
    // label:  

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    val input

    output:   
    path "output_file.txt"

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env bash

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    input = "$input"
    
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