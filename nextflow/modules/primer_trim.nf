process PRIMER_TRIM {
    def module_name = "primer_trim"
    tag "$meta.pcr_primers; $meta.sample_id"
    // label
    // container "pegi3s/cutadapt:latest"
    container "thatdnaguy/cutadapt:v4.7_02"

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_trim_R{1,2}.fastq.gz"),    emit: reads
    path("*_readsout.csv"),                             emit: read_tracking

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
        ${meta.fcid}

    """

}