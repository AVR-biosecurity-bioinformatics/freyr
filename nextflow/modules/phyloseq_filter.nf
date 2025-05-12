process PHYLOSEQ_FILTER {
    def module_name = "phyloseq_filter"
    tag "$pcr_primers"
    label "phyloseq"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(pcr_primers), path(ps), val(loci_params)

    output:
    path("*.csv"),                                                          emit: csvs
    tuple val(pcr_primers), path("ps_filtered_*.rds"), val(loci_params),    emit: ps 
    path("*.fasta"),                                                        emit: asv_fasta
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
    ps =                    "${ps}"
    target_kingdom =        "${loci_params.target_kingdom}"
    target_phylum =         "${loci_params.target_phylum}"
    target_class =          "${loci_params.target_class}"
    target_order =          "${loci_params.target_order}"
    target_family =         "${loci_params.target_family}"
    target_genus =          "${loci_params.target_genus}"
    target_species =        "${loci_params.target_species}"
    min_sample_reads =      "${loci_params.min_sample_reads}"
    min_taxa_reads =        "${loci_params.min_taxa_reads}"
    min_taxa_ra =           "${loci_params.min_taxa_ra}"
    
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
    if ("${params.rdata}" == "true") { save.image(file = "${task.process}.rda") } 
    })

    """

}