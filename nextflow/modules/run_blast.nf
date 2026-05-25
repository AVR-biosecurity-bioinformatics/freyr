process RUN_BLAST {
    def process_name = "run_blast"
    tag "$primers"
    label "medium"
    container "ncbi/blast:2.17.0"

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast_db), val(process_params)

    output:
    tuple val(primers), val(read_group), path(fasta), path("blast.tsv"),                     emit: blast_tsv
    

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    """
    #!/bin/bash

    ### run module code
    bash ${process_name}.sh \
        ${primers} \
        ${read_group} \
        "${fasta}" \
        "${process_params.run_blast}" \
    
    """

}