process SPLIT_LOCI {
    def module_name = "split_loci"
    tag "$sample_primers"
    label "medium"
    container "thatdnaguy/cutadapt:v4.7_02"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)
    val(seq_type)
    val(paired)
    val(primer_error_rate)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("*_R{0,1,2}.fastq.gz"),     emit: reads
    path("*_readsin.csv"),                                                                                  emit: input_counts
    path("*_readsout.csv"),                                                                                 emit: read_tracking

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
        ${process_params.for_primer_seq} \
        ${process_params.rev_primer_seq} \
        ${primers} \
        ${process_params.locus} \
        ${sample_primers} \
        ${read_group} \
        ${primer_error_rate} \
        ${seq_type} \
        ${paired}
    
    """

}