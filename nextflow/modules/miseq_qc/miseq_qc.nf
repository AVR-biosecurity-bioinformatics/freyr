process MISEQ_QC {
    tag "$read_group"
    label "small"
    label "env_r_core"
    publishDir "${params.outdir}/modules/miseq_qc", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${params.outdir}/results/qc/miseq_qc", mode: 'copy'

    input:
    val(read_group) 
    val(miseq_dir)

    output:
    path("*.{pdf,csv}")

    script:
    """

    Rscript --vanilla "${moduleDir}/miseq_qc.R" \
        --process_name "miseq_qc" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --read_group "$read_group" \
        --miseq_dir "$miseq_dir"

    """
}