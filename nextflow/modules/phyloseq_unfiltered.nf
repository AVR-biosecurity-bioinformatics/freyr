process PHYLOSEQ_UNFILTERED {
    def module_name = "phyloseq_unfiltered"
    tag "$pcr_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), path(taxtab), path(seqtab_list), path(samdf_locus), val(loci_params)

    output:
    path("*.csv"),                                                          emit: csvs
    tuple val(pcr_primers), path("ps_unfiltered_*.rds"), val(loci_params),  emit: ps 
    path("*.fasta"),                                                        emit: asv_fasta
    path("accumulation_curve_*.pdf"),                                       emit: acc_curve
    path("*.nwk"),                                                          emit: nwk, optional: true

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
    samdf =                 "${samdf_locus}"
    
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