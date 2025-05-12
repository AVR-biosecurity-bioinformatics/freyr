process TRAIN_IDTAXA {
    def module_name = "train_idtaxa"
    tag "Whole pipeline"
    label "high"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), path(ref_fasta)
    val(fasta_type)

    output:   
    tuple val(pcr_primers), path("*_idtaxa_db.rds")             , emit: model

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${pcr_primers}"
    ref_fasta =             "${ref_fasta}"
    fasta_type =            "${fasta_type}"

    
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