process JOINT_TAX {
    def module_name = "joint_tax"
    tag "$pcr_primers; $fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(fcid), val(loci_params), path(tax, name: "idtaxa*.csv"), path(blast, name: "blast*.csv")
    
    output:
    tuple val(pcr_primers), val(fcid), path("*_joint.csv"), emit: joint

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid =                  "${fcid}"
    pcr_primers =           "${pcr_primers}"
    target_gene =           "${loci_params.target_gene}"
    idtaxa_list =           "${tax}"
    blast_list =            "${blast}"
    
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