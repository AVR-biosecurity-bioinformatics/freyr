process FILTER_QUALPLOTS {
    tag "$sample_primers"
    label "small"
    label "env_r_dada2"
    publishDir { "${params.outdir}/modules/filter_qualplots" }, mode: 'copy', enabled: params.debug_mode.toString().toBoolean()
    publishDir { "${params.outdir}/results/qc/read_filtering" }, mode: 'copy'

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)
    val(file_suffix)

    output:   
    path("qualplots_{pre,post}_*.pdf")                 , emit: plots



    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    Rscript --vanilla "${moduleDir}/filter_qualplots.R" \
        --process_name "filter_qualplots" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --reads_paths "$reads_paths" \
        --locus "$process_params.locus" \
        --seq_type "$process_params.seq_type" \
        --paired "$process_params.paired" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" \
        --file_suffix "$file_suffix"

    """

}