process SPLIT_LOCI {
    def module_name = "split_loci"
    tag "$sample_primers"
    label "medium"
    container "thatdnaguy/cutadapt:v4.7_02"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)

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
        ${process_params.primer_error_rate} \
        ${process_params.seq_type} \
        ${process_params.paired}
    
    """

}