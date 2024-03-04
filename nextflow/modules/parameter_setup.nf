process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    // tag 
    // label

    input:
    path(samdf)
    path(params)

    output:  
    path("input_samdf.rds") ,       emit: input_samdf
    path("params_df.rds") ,         emit: params_df

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    # source functions, themes and load packages
    ### TODO: add these source commands to their own sourced script in "/bin"
    source("${projectDir}/jack_notes/.Rprofile")
    source("${projectDir}/bin/functions.R")
    source("${projectDir}/bin/themes.R")
    source("${projectDir}/bin/_targets_packages.R")
    
    print("${projectDir}")
    print("${samdf}")

    # import inputs as objects
    input_samdf <- read_delim('${samdf}', show_col_types=F, , col_types = cols(.default = "c"))
    input_params <- read_delim('${params}', show_col_types=F, col_types = cols(.default = "c"))

    if(!exists("input_samdf") || !exists("input_params")) {
    stop("Need input_samdf path AND input_params objects to be present. Check module input paths.")
    }
    
    # run module code
    source("${projectDir}/bin/$module_script")

    """

}