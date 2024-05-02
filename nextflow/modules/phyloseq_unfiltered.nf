process PHYLOSEQ_UNFILTERED {
    def module_name = "phyloseq_unfiltered"
    tag "Whole dataset"
    // label:  

    input:
    tuple val(pcr_primers), path(taxtab), path(seqtab_list), val(loci_params)
    path(samdf_original)

    output:
    path("*.csv")
    path("*.rds")
    path("*.fasta")
    path("*.nwk"), optional: true

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    pcr_primers =           "${pcr_primers}"
    taxtab =                "${taxtab}"
    seqtab_list =           "${seqtab_list}"
    loci_params =           "${loci_params}"
    samdf =                 "${samdf_original}"
    
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