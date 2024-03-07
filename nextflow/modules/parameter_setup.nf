process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    // tag: 
    // label: 

    input:
    val data_dir
    // tuple // reads
    // path(samdf)
    // path(params)

    output:  
    path("samdf.csv") ,                 emit: samdf
    path("params.csv") ,                emit: params
    path("params_primer.csv"),          emit: params_primer
    path("params_readfilter.csv"),      emit: params_readfilter
    path("params_dada.csv"),            emit: params_dada
    path("params_asvfilter.csv"),       emit: params_asvfilter
    path("params_database.csv"),        emit: params_database
    path("params_ps.csv"),              emit: params_ps


    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    // def data_dir = "${projectDir}/test_data" // this should defined in workflow script
    """
    #!/usr/bin/env Rscript
    # defining Nextflow environment variables as R variables
    projectDir = "${projectDir}"
    data_dir = "${data_dir}"

    ## source functions, themes and load packages from "bin/process_start.R"
    # this only works with sys.source; "projectDir" doesn't mean anything inside script otherwise
    sys.source("${projectDir}/bin/process_start.R", list(projectDir="${projectDir}"))

    ### run module code
    sys.source(
        "${projectDir}/bin/$module_script", # run script
        envir = .GlobalEnv # this allows import of existing objects like projectDir
        )

    """

}