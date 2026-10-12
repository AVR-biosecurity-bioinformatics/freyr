process MAKE_SEQTAB_PAIRED {
    tag "$primers; $read_group"
    label "long"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/make_seqtab_paired", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(mergers_paths, name: 'mergers*.rds')

    output:
    tuple val(primers), val(read_group), path("*_seqtab_tibble.csv"), path("*_seqs.fasta"),    emit: seqtab
    // path("*_readsout.csv"),                                                                             emit: read_tracking

    script:
    """
    
    make_seqtab_paired.R \
        --process_name "make_seqtab_paired" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" \
        --mergers_paths "$mergers_paths" 

    """

}