process MAKE_SEQTAB_PAIRED {
    tag "$primers; $read_group"
    label "long"
    label "env_r_dada2"
    publishDir "${params.outdir}/modules/make_seqtab_paired", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(mergers_paths, name: 'mergers*.rds')

    output:
    tuple val(primers), val(read_group), path("*_seqtab_tibble.csv"), path("*_seqs.fasta"),    emit: seqtab
    // path("*_readsout.csv"),                                                                             emit: read_tracking

    script:
    """
    
    Rscript --vanilla "${moduleDir}/make_seqtab_paired.R" \
        --process_name "make_seqtab_paired" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" \
        --mergers_paths "$mergers_paths" 

    """

}