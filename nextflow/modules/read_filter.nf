process READ_FILTER {
    def process_name = "read_filter"
    tag "$sample_primers"
    label "medium"
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path(reads), val(process_params)

    output:   
    tuple val(primers), val(read_group), val(sample), val(sample_primers), path("*_filter_R{0,1,2}.fastq.gz"),      emit: reads, optional: true
    path("*_readsout.csv"),                                                                                         emit: read_tracking

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy'

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
        --read_min_length "$process_params.read_min_length" \
        --read_max_length "$process_params.read_max_length" \
        --read_max_ee "$process_params.read_max_ee" \
        --read_trunc_length "$process_params.read_trunc_length" \
        --read_trim_left "$process_params.read_trim_left" \
        --read_trim_right "$process_params.read_trim_right" \
        --read_trim_right "$process_params.read_trim_right" \
        --locus "$process_params.locus" \
        --seq_type "$process_params.seq_type" \
        --paired "$process_params.paired" \
        --sample_primers "$sample_primers" \
        --primers "$primers" \
        --read_group "$read_group" 

    """

}