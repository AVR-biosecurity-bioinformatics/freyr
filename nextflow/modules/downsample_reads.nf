process DOWNSAMPLE_READS {
    def process_name = "downsample_reads"
    tag "$sample"
    label "small"
    container "staphb/seqtk:1.4"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads)
    val(seq_type)
    val(paired)
    val(downsample_reads)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("down_*"), emit: reads

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    bash ${process_name}.sh \
        "${reads_paths}" \
        ${seq_type} \
        ${paired} \
        ${downsample_reads}
    
    """

}