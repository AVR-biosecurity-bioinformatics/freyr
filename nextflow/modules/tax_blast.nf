process TAX_BLAST {
    def module_name = "tax_blast"
    tag "$fcid; $pcr_primers"
    // label:  

    input:
    tuple val(fcid), val(pcr_primers), val(meta), path(seqtab)

    output:
    tuple val(fcid), val(pcr_primers), val(meta), path("*_blast.rds"), emit: blast
    tuple val(fcid), val(pcr_primers), path("*_blast_spp_low.rds"), emit: blast_assignment
    tuple val(fcid), val(pcr_primers), path("n_ranks.txt"), emit: n_ranks

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript

    ### defining Nextflow environment variables as R variables
    ## input channel variables
    fcid =                  "${fcid}"
    pcr_primers =           "${pcr_primers}"
    seqtab =                "${seqtab}"
    target_gene =           "${meta.target_gene}"
    ref_fasta =             "${meta.ref_fasta}"
    blast_min_identity =    "${meta.blast_min_identity}"
    blast_min_coverage =    "${meta.blast_min_coverage}"
    run_blast =             "${meta.run_blast}"


    
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