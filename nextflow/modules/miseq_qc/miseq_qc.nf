process MISEQ_QC {
    tag "$read_group"
    label "small"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/miseq_qc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/miseq_qc", mode: 'copy'

    input:
    val(read_group) 
    val(miseq_dir)

    output:
    path("*.{pdf,csv}")

    script:
    """

    miseq_qc.R \
        --process_name "miseq_qc" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --read_group "$read_group" \
        --miseq_dir "$miseq_dir"

    """
}