process FILTER_SEQTAB {
    def module_name = "filter_seqtab"
    tag "$pcr_primers; $fcid"
    // label:  
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(fcid), val(meta), path(seqtab)

    output:
    tuple val(pcr_primers), val(fcid), val(meta), path("*_seqtab.cleaned.rds"),         emit: seqtab
    tuple val(pcr_primers), val(fcid), val(meta), path("*_ASV_cleanup_summary.csv"),    emit: csv
    tuple val(pcr_primers), val(fcid), val(meta), path("*_ASV_cleanup_summary.pdf"),    emit: plot
    path("*_readsout.csv"),                                                             emit: read_tracking


    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid =              "${fcid}"
    pcr_primers =       "${pcr_primers}"
    seqtab =            "${seqtab}"
    asv_min_length =    "${meta.asv_min_length}"
    asv_max_length =    "${meta.asv_max_length}"
    phmm =              "${meta.phmm}"
    coding =            "${meta.coding}"
    genetic_code =      "${meta.genetic_code}"
    for_primer_seq =    "${meta.for_primer_seq}"
    rev_primer_seq =    "${meta.rev_primer_seq}"

    
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