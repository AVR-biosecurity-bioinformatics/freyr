process DADA_MERGEREADS {
    def module_name = "dada_mergereads"
    tag "$meta.sample_id; $pcr_primers"
    // label:  

    input:
    tuple val(sample_id), val(fcid), val(pcr_primers), val(meta), path(reads), path(seqs)

    output:
    tuple val(sample_id), val(fcid), val(pcr_primers), val(meta), path("*_seqtab.rds"), emit: mergers

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    sample_id =         "${meta.sample_id}"
    fcid =              "${fcid}"
    pcr_primers =       "${pcr_primers}"
    reads_F =           "${reads[0]}"
    reads_R =           "${reads[1]}"
    seqs_F =            "${seqs[0]}"
    seqs_R =            "${seqs[1]}"
    concat_unmerged =   "${meta.concat_unmerged}"
    
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