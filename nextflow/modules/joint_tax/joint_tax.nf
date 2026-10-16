process JOINT_TAX {
    tag "$primers; $read_group"
    label "small"
    label "env_r_core"
    publishDir { "${params.outdir}/modules/joint_tax" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()

    input:
    tuple val(primers), val(read_group), path(tax, name: "idtaxa*.csv"), path(blast, name: "blast*.csv")
    
    output:
    tuple val(primers), val(read_group), path("*_joint.csv"), emit: joint

    script:
    """
    
    Rscript --vanilla "${moduleDir}/joint_tax.R" \
        --process_name "joint_tax" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --primers "$primers" \
        --read_group "$read_group" \
        --idtaxa_list "$tax" \
        --blast_list "$blast" 
        
    """

}