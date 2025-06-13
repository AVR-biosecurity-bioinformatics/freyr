process TAX_BLAST {
    def module_name = "tax_blast"
    tag "$primers; $read_group"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), path(fasta), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_blast.csv"),                     emit: blast
    tuple val(primers), val(read_group), path("*_blast_spp_low.rds"),             emit: blast_assignment
    tuple val(primers), val(read_group), path("*_n_ranks.txt"),                   emit: n_ranks

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
    ref_fasta =             "${process_params.ref_fasta}"
    blast_min_identity =    "${process_params.blast_min_identity}"
    blast_min_coverage =    "${process_params.blast_min_coverage}"
    run_blast =             "${process_params.run_blast}"

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