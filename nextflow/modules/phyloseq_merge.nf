process PHYLOSEQ_MERGE {
    def module_name = "phyloseq_merge"
    tag "Whole dataset"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(ps_unfiltered)
    path(ps_filtered)

    output:
    path("*.csv")
    path("*.rds")
    path("*.fasta")
    path("*.nwk"),          optional: true
    path("*_readsout.csv"), emit: read_tracking

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    ps_unfiltered =                  "${ps_unfiltered}"
    ps_filtered =                    "${ps_filtered}"
    
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