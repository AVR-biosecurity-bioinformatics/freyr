process TEMPLATE {
    def module_name = "template"

    // tag 
    // label 

    input:

    output:   

    publishDir "${projectDir}/output/modules/${module_name}"

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    # source functions, themes and load packages
    source("${projectDir}/bin/functions.R")
    source("${projectDir}/bin/themes.R")
    source("${projectDir}/bin/_targets_packages.R")
    
    # import inputs as objects
    samdf_file <- read_delim($samdf_file, show_col_types=F)
    params_file <- read_delim($loci_params, show_col_types=F)

    if(!exists("samdf_file") || !exists("params_file")) {
    stop("Need samdf path AND params_file objects to be present. Check module input paths.")

    # run module code
    source("${projectDir}/bin/$module_script")

    """

}