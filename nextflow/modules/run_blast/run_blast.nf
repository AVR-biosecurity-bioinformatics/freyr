process RUN_BLAST {
    tag "$primers"
    label "medium"
    label "env_blast"
    publishDir { "${params.outdir}/modules/run_blast" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast_db), val(process_params)

    output:
    tuple val(primers), val(read_group), path(fasta), path("blast.tsv"),                     emit: blast_tsv

    script:
    """
    #!/bin/bash

    ### run module code
    ${moduleDir}/run_blast.sh \
        ${primers} \
        ${read_group} \
        "${fasta}" \
        "${process_params.run_blast}" \
    
    """

}