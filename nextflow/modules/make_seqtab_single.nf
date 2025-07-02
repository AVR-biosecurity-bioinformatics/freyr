process MAKE_SEQTAB_SINGLE {
    def process_name = "make_seqtab_single"
    tag "$primers; $read_group"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path(seqs), val(concat_unmerged)

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
        --reads "$reads" \
        --seqs "$seqs" \
        --concat_unmerged "$concat_unmerged"

    """

}