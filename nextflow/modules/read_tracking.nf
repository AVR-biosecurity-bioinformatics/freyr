process READ_TRACKING {
    def module_name = "read_tracking"
    tag "Whole dataset"
    // label:  
    // cache false 

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