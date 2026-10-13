process NANOPLOT {
    tag "$sample"
    label "high"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/nanoplot", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/nanoplot", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*report.html"),     emit: report

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/nanoplot.sh \
        "${reads_paths}" \
        ${sample} \
        ${seq_type} \
        ${paired}
    
    """

}