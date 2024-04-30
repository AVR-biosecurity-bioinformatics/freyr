process TAX_SUMMARY {
    def module_name = "tax_summary"
    tag "$fcid; $pcr_primers"
    // label:  

    input:
    tuple val(pcr_primers), val(fcid), path(tax), path(ids), val(target_gene), path(joint), path(idtaxa_db), path(ref_fasta)

    output:
    tuple val(pcr_primers), val(fcid), val(target_gene), path("*_taxonomic_assignment_summary.rds"), emit: rds
    path("*_taxonomic_assignment_summary.csv"), emit: csv

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
    tax =                   "${tax}"
    ids =                   "${ids}"
    target_gene =           "${target_gene}"
    joint_file =            "${joint}"
    idtaxa_db =             "${idtaxa_db}"
    ref_fasta =             "${ref_fasta}"
    
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