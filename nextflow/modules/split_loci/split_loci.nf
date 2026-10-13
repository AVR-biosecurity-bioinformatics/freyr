process SPLIT_LOCI {
    tag "$sample_primers"
    label "medium"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/split_loci", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("*_R{0,1,2}.fastq.gz"),     emit: reads
    path("*_readsin.csv"),                                                                                  emit: input_counts
    path("*_readsout.csv"),                                                                                 emit: read_tracking

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/split_loci.sh \
        "${reads_paths}" \
        ${process_params.for_primer_seq} \
        ${process_params.rev_primer_seq} \
        ${primers} \
        ${process_params.locus} \
        ${sample_primers} \
        ${read_group} \
        ${process_params.primer_error_rate} \
        ${process_params.seq_type} \
        ${process_params.paired}
    
    """

}