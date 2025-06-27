process ASSIGNMENT_PLOT {
    def module_name = "assignment_plot"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast, name: "blast*.rds"), path(tax), val(process_params)

    output:
    path("*_taxonomic_assignment_summary.pdf"),                 emit: plot
    tuple val(primers), val(read_group), path("*_joint.rds"),   emit: joint

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    read_group =            "${read_group}"
    primers =               "${primers}"
    fasta =                 "${fasta}"
    blast_list =            "${blast}"
    tax =                   "${tax}"
    idtaxa_db =             "${process_params.idtaxa_db}"
    ref_fasta =             "${process_params.ref_fasta}"
    
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