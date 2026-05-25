process MAKE_BLAST_DATABASE {
    def process_name = "make_blast_database"
    tag "$primers"
    label "medium"
    container "ncbi/blast:2.17.0"

    input:
    tuple val(primers), val(process_params)

    output:
    tuple val(primers), path("*.blast_db.*"),                     emit: blast_db
    

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    """
    #!/bin/bash

    ### run module code
    bash ${process_name}.sh \
        ${primers} \
        "${process_params.ref_fasta}"
    
    """

}