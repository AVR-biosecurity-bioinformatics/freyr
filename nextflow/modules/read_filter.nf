process READ_FILTER {
    def module_name = "read_filter"
    tag "$meta.sample_id; $meta.target_gene"
    // label:  

    input:
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("R{1,2}_out.fastq.gz"),           emit: reads


    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fwd_reads = "${reads[0]}"
    rev_reads = "${reads[1]}"
    read_min_length = "${meta.read_min_length}"
    read_max_length = "${meta.read_max_length}"
    read_max_ee = "${read_max_ee}"
    read_trunc_length = "${read_trunc_length}"
    read_trim_left = "${read_trim_left}"
    read_trim_right = "${read_trim_right}"
    
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