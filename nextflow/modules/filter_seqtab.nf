process FILTER_SEQTAB {
    def module_name = "filter_seqtab"
    tag "$fcid; $pcr_primers"
    // label:  

    input:
    tuple val(fcid), val(pcr_primers), val(meta), path(seqtab)

    output:

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
    asv_max_length =    "${asv_max_length}"
    phmm =              "${meta.phmm}"
    coding =            "${meta.coding}"
    genetic_code =      "${genetic_code}"
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

    """

}