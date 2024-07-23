process SPLIT_LOCI {
    def module_name = "split_loci"
    tag "$meta.pcr_primers; $meta.sample_id"
    label "medium"
    // container "nanozoo/bbmap:38.86--9ebcbfa"
    container "thatdnaguy/cutadapt:v4.7_02"

    input:
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_R{1,2}.fastq.gz"),     emit: reads
    // path("split_loci_*.txt")
    path("*_readsin.csv"),                          emit: input_counts
    path("*_readsout.csv"),                         emit: read_tracking

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/bin/bash

    ### run module code
    bash ${module_name}.sh \
        ${reads[0]} \
        ${reads[1]} \
        ${meta.for_primer_seq} \
        ${meta.rev_primer_seq} \
        ${meta.pcr_primers} \
        ${meta.target_gene} \
        ${meta.sample_id} \
        ${meta.fcid} \
        ${params.primer_error_rate}
    
    """

}