process FILTER_QUALPLOTS {
    tag "$sample_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"
    publishDir "${launchDir}/output/modules/filter_qualplots", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/qc/read_filtering", mode: 'copy'

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)
    val(file_suffix)

    output:   
    path("qualplots_{pre,post}_*.pdf")                 , emit: plots



    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """
    filter_qualplots.R \
        --process_name "filter_qualplots" \
        --projectDir "$projectDir" \
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