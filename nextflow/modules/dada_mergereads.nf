process DADA_MERGEREADS {
    def module_name = "dada_mergereads"
    // tag "$fcid; $pcr_primers"
    // label:  

    input:
    tuple val(sample_id), val(fcid), val(pcr_primers), val(meta), path(seqF), path(seqR)

    output:

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
    seqF =              "${seqF}"
    seqR =              "${seqR}"
    
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