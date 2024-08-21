process DOWNSAMPLE_READS {
    def module_name = "downsample_reads"
    tag "$meta.pcr_primers; $meta.sample_id"
    label "medium"
    container "staphb/seqtk:1.4"

    input:
    tuple val(meta), path(reads)
    val(seq_type)
    val(paired)
    val(downsample_reads)

    output:   
    tuple val(meta), path("down*")            , emit: reads

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
        ${seq_type} \
        ${paired} \
        ${downsample_reads}
    
    """

}