process RUN_BLAST {
    tag "$primers"
    label "medium"
    container "ncbi/blast:2.17.0"
    publishDir "${launchDir}/output/modules/run_blast", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast_db), val(process_params)

    output:
    tuple val(primers), val(read_group), path(fasta), path("blast.tsv"),                     emit: blast_tsv

    script:
    """
    #!/bin/bash

    ### run module code
    bash run_blast.sh \
        ${primers} \
        ${read_group} \
        "${fasta}" \
        "${process_params.run_blast}" \
    
    """

}