process FILTER_QUALPLOTS {
    def module_name = "filter_qualplots"
    tag "$meta.fcid; $meta.sample_id"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(meta), path(reads)
    val(seq_type)
    val(paired)
    val(file_suffix)

    output:   
    path("*_qualplots.pdf")                 , emit: plots, optional: true

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"

    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    reads_paths =       "${reads_paths}"
    sample_id =         "${meta.sample_id}"
    fcid =              "${meta.fcid}"
    target_gene =       "${meta.target_gene}"
    pcr_primers =       "${meta.pcr_primers}"
    seq_type =          "${seq_type}"
    paired =            "${paired}"
    file_suffix =       "${file_suffix}"
    
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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}_${task.index}.rda") } 
    })
    
    """

}