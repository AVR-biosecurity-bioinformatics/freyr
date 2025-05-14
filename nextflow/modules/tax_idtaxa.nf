process TAX_IDTAXA {
    def module_name = "tax_idtaxa"
    tag "$pcr_primers; $fcid"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(fcid), val(loci_params), path(fasta)
    //tuple val(pcr_primers2), val(fcid2), val(loci_params2), path(fasta)

    output:
    tuple val(pcr_primers), val(fcid), val(loci_params), path("*_idtaxa_tax.rds"), emit: tax
    tuple val(pcr_primers), val(fcid), val(loci_params), path("*_idtaxa_ids.rds"), emit: ids

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
    idtaxa_confidence = "${loci_params.idtaxa_confidence}"
    idtaxa_db =         "${loci_params.idtaxa_db}"
    fasta =             "${fasta}"

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