process READ_TRACKING {
    def process_name = "read_tracking"
    tag "Whole dataset"
    label "small"
    // cache false 
    container "jackscanlan/piperline-multi:0.0.1"

    input:
    path(rt_samples)
    path(rt_group)
    path(samplesheet_split)

    output:
    path("*.csv")
    path("read_tracker.csv"),           emit: csv
    path("read_tracker.pdf"),           emit: plot
    path("*.pdf")

    publishDir "${launchDir}/output/modules/${process_name}", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"

    publishDir "${launchDir}/output/results/read_tracking", mode: 'copy'


    // when: 

    script:
    """
    
    ${process_name}.R \
        --process_name "$process_name" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --rt_samples "$rt_samples" \
        --rt_group "$rt_group" \
        --samplesheet_split_file "$samplesheet_split" 
        
    """

}