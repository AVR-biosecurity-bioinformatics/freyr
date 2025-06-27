process DOWNSAMPLE_READS {
    def module_name = "downsample_reads"
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