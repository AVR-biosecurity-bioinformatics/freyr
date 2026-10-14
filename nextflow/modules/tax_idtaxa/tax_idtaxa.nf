process TAX_IDTAXA {
    tag "$primers; $read_group"
    label "medium"
    label "env_r_phyloseq"
    publishDir "${params.outdir}/modules/tax_idtaxa", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    input:
    tuple val(primers), val(read_group), path(fasta), val(process_params)

    output:
    tuple val(primers), val(read_group), path("*_idtaxa_tax.csv"), emit: tax
    tuple val(primers), val(read_group), path("*_idtaxa_ids.rds"), emit: ids

    script:
    """
    
    Rscript --vanilla "${moduleDir}/tax_idtaxa.R" \
        --process_name "tax_idtaxa" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --idtaxa_confidence "$process_params.idtaxa_confidence" \
        --idtaxa_db "$process_params.idtaxa_db" \
        --fasta "$fasta" 
        
    """

}