process TAX_SUMMARY {
    def module_name = "tax_summary"
    tag "$pcr_primers; $fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(fcid), path(tax), path(ids), path(joint), val(loci_params)

    output:
    tuple val(pcr_primers), val(fcid), val(loci_params), path("*_taxonomic_assignment_summary.rds"), emit: rds
    path("*_taxonomic_assignment_summary.csv"), emit: csv

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${pcr_primers}"
    fcid =                  "${fcid}"   
    tax =                   "${tax}"
    ids =                   "${ids}"
    joint_file =            "${joint}"
    target_gene =           "${loci_params.target_gene}"
    idtaxa_db =             "${loci_params.idtaxa_db}"
    ref_fasta =             "${loci_params.ref_fasta}"
    
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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}_${task.index}.rda") } 
    })

    """

}