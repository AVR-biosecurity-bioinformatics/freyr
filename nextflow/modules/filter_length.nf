process FILTER_LENGTH {
    def module_name = "filter_length"
    tag "$pcr_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(meta), path(seqtab_tibble_list), path(fasta_list)

    output:
    tuple val(pcr_primers), path("*_length_filter.csv"),            emit: tibble

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
        
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${pcr_primers}"
    seqtab_tibble_list =    "${seqtab_tibble_list}"
    fasta_list =            "${fasta_list}"
    asv_min_length =        "${meta.asv_min_length}"
    asv_max_length =        "${meta.asv_max_length}"
    
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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """

}