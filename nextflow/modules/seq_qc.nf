process SEQ_QC {
    def module_name = "seq_qc"
    // tag:
    // label: 

    input:
    path sample_info
    path loci_params

    output:

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'


    // when:

    script:
    def module_script = "${module_name}.R"
    """
    #!/usr/bin/env Rscript
    
    # source functions, themes and load packages from "bin/process_start.R"
    # this only works this way; "projectDir" doesn't mean anything inside script
    sys.source("${projectDir}/bin/process_start.R", list(projectDir="${projectDir}"))

    # run module code
    source("${projectDir}/bin/$module_script")

    """
}