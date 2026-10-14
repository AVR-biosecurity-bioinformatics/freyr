process FASTQC {
    tag "$sample"
    label "small"
    label "env_fastqc"
    publishDir { "${params.outdir}/modules/fastqc" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/qc/fastqc" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*.html"),     emit: report

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/fastqc.sh \
        "${reads_paths}" \
        ${task.memory.mega} \
        ${seq_type} \
        ${paired}
    
    """

}