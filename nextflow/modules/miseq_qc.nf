process MISEQ_QC {
    def process_name = "miseq_qc"
    tag "$read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    val(read_group) 
    val(miseq_dir)

    output:
    path("*.{pdf,csv}")

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/qc/miseq_qc", mode: 'copy'

    // when:

    script:
    """

    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --read_group "$read_group" \
        --miseq_dir "$miseq_dir"

    """
}