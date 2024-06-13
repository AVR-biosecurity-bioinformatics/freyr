process READ_FILTER {
    def module_name = "read_filter"
    tag "$meta.pcr_primers; $meta.sample_id"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_filter_R{1,2}.fastq.gz"),           emit: reads, optional: true
    path("*_readsout.csv"),                                      emit: read_tracking

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
    fcid =              "${meta.fcid}"
    
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