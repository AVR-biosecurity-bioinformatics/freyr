process PRIMER_TRIM {
    def module_name = "primer_trim"
    tag "$meta.pcr_primers; $meta.sample_id"
    // label:  

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_trim_R{1,2}.fastq.gz"), emit: reads
    tuple val(module_name), val(meta.pcr_primers), val(meta.fcid), val(meta.sample_id), path("R1_output.txt"), path("R2_output.txt"), emit: output_counts

    publishDir "${projectDir}/output/modules/${module_name}", mode: 'copy'

    // when: 

    script:
    def module_script = "${module_name}.sh"
    """
    #!/usr/bin/bash

    ### run module code
    bash ${module_name}.sh \
        ${reads[0]} \
        ${reads[1]} \
        ${meta.for_primer_seq} \
        ${meta.rev_primer_seq} \
        ${meta.pcr_primers} \
        ${meta.target_gene} \
        ${meta.sample_id} 

    """

}