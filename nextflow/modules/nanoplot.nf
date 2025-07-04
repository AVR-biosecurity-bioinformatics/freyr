process NANOPLOT {
    def process_name = "nanoplot"
    tag "$sample"
    label "high"
    container "nanozoo/nanoplot:1.42.0--547049c"

    input:
    tuple val(sample), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*report.html"),     emit: report

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/qc/nanoplot", mode: 'copy'

    // when: 

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    bash ${process_name}.sh \
        "${reads_paths}" \
        ${sample} \
        ${seq_type} \
        ${paired}
    
    """

}