process NANOPLOT {
    def module_name = "nanoplot"
    tag "$meta.pcr_primers; $meta.sample_id"
    label "high"
    container "nanozoo/nanoplot:1.42.0--547049c"

    input:
    tuple val(meta), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*report.html"),     emit: report

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"

    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    bash ${module_name}.sh \
        "${reads_paths}" \
        ${meta.sample_id} \
        ${seq_type} \
        ${paired}
    
    """

}