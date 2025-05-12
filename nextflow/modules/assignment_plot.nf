process ASSIGNMENT_PLOT {
    def module_name = "assignment_plot"
    tag "$pcr_primers; $fcid"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), val(fcid), val(loci_params), path(seqtab), path(blast), path(tax)

    output:
    path("*_taxonomic_assignment_summary.pdf"), emit: plot
    tuple val(pcr_primers), val(fcid), val(loci_params), path("*_joint.rds"), emit: joint

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
    seqtab =                "${seqtab}"
    blast =                 "${blast}"
    tax =                   "${tax}"
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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """

}