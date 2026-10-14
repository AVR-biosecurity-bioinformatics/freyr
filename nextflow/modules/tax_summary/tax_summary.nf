process TAX_SUMMARY {
    tag "$primers; $read_group"
    label "small"
    label "env_r_core"
    publishDir { "${params.outdir}/modules/tax_summary" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/taxonomic_assignment" }, pattern: 'taxonomic_assignment_summary_*.csv', mode: 'copy'

    input:
    tuple val(primers), val(read_group), path(fasta), path(tax, name: "idtaxa*.csv"), path(ids, name: "idtaxa*.rds"), path(joint), val(process_params)

    output:
    tuple val(primers), val(read_group), path("taxonomic_assignment_summary_*.rds"), emit: rds
    tuple val(primers), val(read_group), path("taxonomic_assignment_summary_*.csv"), emit: csv

    script:
    """
    
    Rscript --vanilla "${moduleDir}/tax_summary.R" \
        --process_name "tax_summary" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --fasta "$fasta" \
        --tax_list "$tax" \
        --ids_list "$ids" \
        --joint_file "$joint" \
        --locus "$process_params.locus" \
        --idtaxa_db "$process_params.idtaxa_db" \
        --ref_fasta "$process_params.ref_fasta"
        
    """

}