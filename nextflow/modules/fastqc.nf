process FASTQC {
    def process_name = "fastqc"
    tag "$sample"
    label "medium"
    container "staphb/fastqc:0.12.1"

    input:
    tuple val(sample), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*.html"),     emit: report

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/qc/fastqc", mode: 'copy'


    // when: 

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    bash ${process_name}.sh \
        "${reads_paths}" \
        ${task.memory.mega} \
        ${seq_type} \
        ${paired}
    
    """

}