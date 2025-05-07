process FASTQC {
    def module_name = "fastqc"
    tag "$sample_id"
    label "medium"
    container "staphb/fastqc:0.12.1"

    input:
    tuple val(sample_id), path(reads)
    val(seq_type)
    val(paired)

    output:   
    path("*.html"),     emit: report

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
        ${task.memory.mega} \
        ${seq_type} \
        ${paired}
    
    """

}