process FILTER_QUALPLOTS {
    def process_name = "filter_qualplots"
    tag "$sample_primers"
    label "small"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)
    val(file_suffix)

    output:   
    path("*_qualplots.pdf")                 , emit: plots, optional: true

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    // when: 

    script:
    def reads_paths = reads.join(";") // concatenate reads into a single string
    """

    ${process_name}.R \
        --process_name "$process_name" \
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