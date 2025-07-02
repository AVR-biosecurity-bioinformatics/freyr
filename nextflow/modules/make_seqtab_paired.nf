process MAKE_SEQTAB_PAIRED {
    def process_name = "make_seqtab_paired"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(readsF), path(readsR), path(seqsF), path(seqsR), val(concat_unmerged)

    output:
    tuple val(primers), val(read_group), path("*_seqtab_tibble.csv"), path("*_seqs.fasta"),    emit: seqtab
    path("*_readsout.csv"),                                                                             emit: read_tracking

    publishDir "${projectDir}/output/modules/${process_name}", mode: 'copy'

    // when: 

    script:
    """

    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" \
        --reads_F "$readsF" \
        --reads_R "$readsR" \
        --seqs_F "$seqsF" \
        --seqs_R "$seqsR" \
        --concat_unmerged "$concat_unmerged"

    """

}