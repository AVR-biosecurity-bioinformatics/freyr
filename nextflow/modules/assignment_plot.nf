process ASSIGNMENT_PLOT {
    def module_name = "assignment_plot"
    tag "$pcr_primers; $fcid"
    // label:  

    input:
    tuple val(pcr_primers), val(fcid), path(seqtab), path(blast), path(tax), val(target_gene), path(idtaxa_db), path(ref_fasta)

    output:
    path("*_taxonomic_assignment_summary.pdf"), emit: plot
    tuple val(pcr_primers), val(fcid), val(target_gene), path("*_joint.rds"), path(idtaxa_db), path(ref_fasta), emit: joint

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
    blast =                 "${blast}"
    tax =                   "${tax}"
    target_gene =           "${target_gene}"
    idtaxa_db =             "${idtaxa_db}"
    ref_fasta =             "${ref_fasta}"
    
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