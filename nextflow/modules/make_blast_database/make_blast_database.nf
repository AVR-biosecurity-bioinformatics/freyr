process MAKE_BLAST_DATABASE {
    tag "$primers"
    label "medium"
    label "env_blast"
    publishDir { "${params.outdir}/modules/make_blast_database" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), val(process_params)

    output:
    tuple val(primers), path("*.blast_db.*"),                     emit: blast_db

    script:
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/make_blast_database.sh \
        ${primers} \
        "${process_params.ref_fasta}"
    
    """

}