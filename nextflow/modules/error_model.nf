process ERROR_MODEL {
    def module_name = "error_model"
    // tag: 
    // label:  

    input:
    val direction
    tuple val(meta), path(reads)


    output:   
    path "output_file.txt"

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fwd_reads =         "${reads[0]}"
    rev_reads =         "${reads[1]}"
    read_min_length =   "${meta.read_min_length}"
    read_max_length =   "${meta.read_max_length}"
    read_max_ee =       "${meta.read_max_ee}"
    read_trunc_length = "${meta.read_trunc_length}"
    read_trim_left =    "${meta.read_trim_left}"
    read_trim_right =   "${meta.read_trim_right}"
    sample_id =         "${meta.sample_id}"
    target_gene =       "${meta.target_gene}"
    pcr_primers =       "${meta.pcr_primers}"
    direction =         "${direction}"
    
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