process SPLIT_LOCI {
    def module_name = "split_loci"
    tag "$meta.sample_id; $meta.target_gene"
    // label:  

    input:
    // input read pairs, primer seqs and locus name (to use in output file names)
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("${meta.sample_id}_${meta.target_gene}_${meta.pcr_primers}_R*.fastq.gz"), emit: reads
    path("split_loci_*.txt")

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