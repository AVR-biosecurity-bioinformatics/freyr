process READ_TRACKING {
    tag "Whole dataset"
    label "small"
    // cache false 
    conda "${moduleDir}/environment.yml"
    publishDir "${launchDir}/output/modules/read_tracking", mode: 'copy', enabled: "${ params.debug_mode ? true : false }"
    publishDir "${launchDir}/output/results/read_tracking", mode: 'copy'
    
    input:
    path(rt_samples, name: 'sample*.csv')
    path(rt_group, name: 'group*.csv')
    path(samplesheet_split)

    output:
    path("sample_tibble.csv")
    path("group_tibble.csv")
    path("read_tracker.csv"),           emit: csv
    path("read_tracker.pdf"),           emit: plot
    path("*.pdf")

    script:
    """
    
    read_tracking.R \
        --process_name "read_tracking" \
        --projectDir "$projectDir" \
        --cpus "$task.cpus" \
        --rdata "$params.rdata" \
        --rt_samples "$rt_samples" \
        --rt_group "$rt_group" \
        --samplesheet_split_file "$samplesheet_split" 
        
    """

}