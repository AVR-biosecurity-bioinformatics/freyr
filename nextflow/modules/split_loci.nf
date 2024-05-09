process SPLIT_LOCI {
    def module_name = "split_loci"
    tag "$meta.pcr_primers; $meta.sample_id"
    // label:  

    input:
    tuple val(meta), path(reads)

    output:   
    tuple val(meta), path("*_R{1,2}.fastq.gz"), emit: reads
    path("split_loci_*.txt")
    tuple val("input"), val(meta.pcr_primers), val(meta.fcid), val(meta.sample_id), path("readsin.csv"), emit: input_counts
    path("readsout.csv"), emit: output_counts

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
        ${meta.sample_id} \
        ${meta.fcid}
    
    """

}