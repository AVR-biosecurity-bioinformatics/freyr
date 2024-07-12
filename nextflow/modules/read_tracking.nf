process READ_TRACKING {
    def module_name = "read_tracking"
    tag "Whole dataset"
    label "small"
    // cache false 
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(rt_samples)
    path(rt_group)

    output:
    path("*.csv")
    path("read_tracker.csv"),           emit: csv
    path("read_tracker.pdf"),           emit: plot
    path("*.pdf")

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    rt_samples =         "${rt_samples}"
    rt_group =           "${rt_group}"
    
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