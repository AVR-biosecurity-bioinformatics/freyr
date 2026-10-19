process MAKE_SEQTAB_SINGLE {
    tag "$primers; $read_group"
    label "small"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/make_seqtab_single" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), path(seqs), val(concat_unmerged)

    output:
    tuple val(primers), val(read_group), path("*_seqtab_tibble.csv"), path("*_seqs.fasta"),    emit: seqtab
    path("*_readsout.csv"),                                                                             emit: read_tracking

    script:
    """

    Rscript --vanilla "${moduleDir}/make_seqtab_single.R" \
        --process_name="make_seqtab_single" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --sample_primers="$sample_primers" \
        --primers="$primers" \
        --read_group="$read_group" \
        --reads="$reads" \
        --seqs="$seqs" \
        --concat_unmerged="$concat_unmerged"

    """

}