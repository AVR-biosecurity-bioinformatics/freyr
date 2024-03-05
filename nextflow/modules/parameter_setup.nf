process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    // tag: 
    // label: 

    input:
    path(samdf)
    path(params)

    output:  
    path("input_samdf.csv") ,           emit: input_samdf
    path("input_params.csv") ,          emit: input_params
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
    """
    #!/usr/bin/env Rscript

    # source functions, themes and load packages from "bin/process_start.R"
    # this only works this way; "projectDir" doesn't mean anything inside script
    sys.source("${projectDir}/bin/process_start.R", list(projectDir="${projectDir}"))
    
    print("${projectDir}")
    print("${samdf}")

    # import inputs as objects
    input_samdf <- read_delim('${samdf}', show_col_types=F, , col_types = cols(.default = "c"))
    input_params <- read_delim('${params}', show_col_types=F, col_types = cols(.default = "c"))

    if(!exists("input_samdf") || !exists("input_params")) {
    stop("Need input_samdf path AND input_params objects to be present. Check module input paths.")
    }
    
    # create directories required by pipeline if they don't exist
    # (formerly 'create_dirs' target)
    step_validate_folders("${projectDir}")

    # run module code
    source("${projectDir}/bin/$module_script")

    """

}