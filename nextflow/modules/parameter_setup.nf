process PARAMETER_SETUP {
    def module_name = "parameter_setup"
    // tag 
    // label
    stageInMode = 'copy' 

    input:
    path(samdf)
    path(params)

    output:  
    path("input_samdf.csv") ,       emit: input_samdf
    path("params_df.csv") ,         emit: params_df

    publishDir "${projectDir}/output/modules/${module_name}"

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
    samdf_file <- read_delim(${samdf}, show_col_types=F)
    params_file <- read_delim(${params}, show_col_types=F)

    if(!exists("samdf_file") || !exists("params_file")) {
    stop("Need samdf path AND params_file objects to be present. Check module input paths.")

    # run module code
    source("${projectDir}/bin/$module_script")

    """

}