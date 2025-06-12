process FILTER_FRAME {
    def module_name = "filter_frame"
    tag "$primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), path(seqtab_tibble_list), path(fasta_list), val(process_params)

    output:
    tuple val(primers), path("*_frame_filter.csv"),       emit: tibble

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
        
    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${primers}"
    seqtab_tibble_list =    "${seqtab_tibble_list}"
    fasta_list =            "${fasta_list}"
    coding =                "${process_params.coding}"
    genetic_code =          "${process_params.genetic_code}"
    
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