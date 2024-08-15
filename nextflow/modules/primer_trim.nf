process PRIMER_TRIM {
    def module_name = "primer_trim"
    tag "$meta.pcr_primers; $meta.sample_id"
    label "medium"
    container "thatdnaguy/cutadapt:v4.7_02"

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_trim_R{0,1,2}.fastq.gz"),    emit: reads
    path("*_readsout.csv"),                             emit: read_tracking

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"

    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    #!/bin/bash

    ### run module code
    bash ${module_name}.sh \
        ${reads_paths} \
        ${meta.for_primer_seq} \
        ${meta.rev_primer_seq} \
        ${meta.pcr_primers} \
        ${meta.target_gene} \
        ${meta.sample_id} \
        ${meta.fcid} \
        ${params.primer_n_trim} \
        ${params.primer_error_rate}

    """

}