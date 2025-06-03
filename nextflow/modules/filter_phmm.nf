process FILTER_PHMM {
    def module_name = "filter_phmm"
    tag "$pcr_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(meta), path(seqtab_tibble_list), path(fasta_list)

    output:
    tuple val(pcr_primers), path("*_phmm_filter.csv"),                     emit: tibble
    path("*_readsout.csv"),                                                   emit: read_tracking


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
    phmm =                  "${meta.phmm}"
    for_primer_seq =        "${meta.for_primer_seq}"
    rev_primer_seq =        "${meta.rev_primer_seq}"
    coding =                "${meta.coding}"
    
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