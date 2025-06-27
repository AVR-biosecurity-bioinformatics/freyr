process TAX_SUMMARY {
    def module_name = "tax_summary"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), path(tax, name: "idtaxa*.csv"), path(ids, name: "idtaxa*.rds"), path(joint), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_taxonomic_assignment_summary.rds"), emit: rds
    tuple val(primers), val(read_group), path("*_taxonomic_assignment_summary.csv"), emit: csv

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    primers =               "${primers}"
    read_group =            "${read_group}"   
    fasta =                 "${fasta}"
    tax_list =              "${tax}"
    ids_list =              "${ids}"
    joint_file =            "${joint}"
    locus =                 "${process_params.locus}"
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