process MERGE_PAIRED {
    tag "$sample_primers"
    label "medium"
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/merge_paired", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(readsF), path(readsR), path(seqsF), path(seqsR), val(concat_unmerged)

    output:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("mergers.rds"),         emit: mergers
    path("*_readsout.csv"),                                                                             emit: read_tracking

    script:
    """
     
    merge_paired.R \
        --process_name "merge_paired" \
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