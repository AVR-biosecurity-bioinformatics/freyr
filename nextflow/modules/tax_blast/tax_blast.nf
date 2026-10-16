process TAX_BLAST {
    tag "$primers; $read_group"
    label "medium"
    label "env_r_core"
    publishDir { "${params.outdir}/modules/tax_blast" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), val(read_group), path(fasta), path(blast_tsv), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_blast.csv"),                     emit: blast
    tuple val(primers), val(read_group), path("*_blast_spp_low.rds"),             emit: blast_assignment
    // tuple val(primers), val(read_group), path("*_n_ranks.txt"),                   emit: n_ranks

    script:
    """
    
    Rscript --vanilla "${moduleDir}/tax_blast.R" \
        --process_name "tax_blast" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --fasta "$fasta" \
        --blast_tsv "$blast_tsv" \
        --blast_min_identity "$process_params.blast_min_identity" \
        --blast_min_coverage "$process_params.blast_min_coverage" \
        --run_blast "$process_params.run_blast"
        
    """

}