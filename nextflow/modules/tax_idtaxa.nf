process TAX_IDTAXA {
    def module_name = "tax_idtaxa"
    tag "$pcr_primers; $fcid"
    // label:  

    input:
    tuple val(pcr_primers), val(fcid), val(meta), path(seqtab)

    output:
    tuple val(pcr_primers), val(fcid), val(meta), path("*_idtaxa_tax.rds"), emit: tax
    tuple val(pcr_primers), val(fcid), val(meta), path("*_idtaxa_ids.rds"), emit: ids



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
    idtaxa_confidence = "${meta.idtaxa_confidence}"
    idtaxa_db =         "${meta.idtaxa_db}"


    
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