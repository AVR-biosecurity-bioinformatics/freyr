process DADA_MERGEREADS {
    def module_name = "dada_mergereads"
    tag "$pcr_primers; $fcid"
    // label:  

    input:
    tuple val(pcr_primers), val(fcid), val(concat_unmerged), val(meta), path(readsF), path(readsR), path(seqsF), path(seqsR)

    output:
    tuple val(pcr_primers), val(fcid), val(meta), path("*_seqtab.rds"), emit: seqtab
    path("*_readsout.csv"),                                             emit: read_tracking

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
    sample_id =         "${meta.sample_id}"
    reads_F =           "${readsF}"
    reads_R =           "${readsR}"
    seqs_F =            "${seqsF}"
    seqs_R =            "${seqsR}"
    concat_unmerged =   "${concat_unmerged}"
    
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