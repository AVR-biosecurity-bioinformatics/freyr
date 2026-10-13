process DOWNSAMPLE_READS {
    tag "$sample"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/modules/downsample_reads", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads)
    val(seq_type)
    val(paired)
    val(downsample_reads)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("down_*"), emit: reads

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/downsample_reads.sh \
        "${reads_paths}" \
        ${seq_type} \
        ${paired} \
        ${downsample_reads}
    
    """

}