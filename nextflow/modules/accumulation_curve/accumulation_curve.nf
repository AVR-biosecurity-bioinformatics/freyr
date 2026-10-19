process ACCUMULATION_CURVE {
    tag "$primers"
    label "high"
    label "env_r_phyloseq"
    publishDir { "${params.outdir}/modules/accumulation_curve" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/main_outputs/unfiltered" }, mode: 'copy'

    input:
    tuple val(primers), path(seqtab), path(samdf), val(process_params)

    output:
    path("accumulation_curve_*.pdf"),                                          emit: pdf

    script:
    """
    Rscript --vanilla "${moduleDir}/accumulation_curve.R" \
        --process_name="accumulation_curve" \
        --cpus="$task.cpus" \
        --rdata="$params.rdata" \
        --primers="$primers" \
        --seqtab_file="$seqtab" \
        --samdf_file="$samdf" \
        --min_sample_reads="$process_params.min_sample_reads"
        
    """

}