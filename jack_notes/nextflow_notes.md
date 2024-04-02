> Add custom `.Rprofile` to GitHub repo (in `/jack_notes` dir?) that can be sourced at run time. This contains the location of the library that should have all the packages needed for analysis. 


Need to source copy of `.Rprofile` that contains new default library location and option to use `pak` in `renv`. 

    # load modules
    module load Nextflow
    module load R/4.2.0-foss-2021b
    module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
    module load GDAL/3.3.0-foss-2021a
    module load BLAST+/2.11.0-gompi-2020a
    module load Pandoc/2.5
    module load ZeroMQ/4.3.2-GCCcore-9.3.0

    # load module oneline
    module load Nextflow R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0 
    
    ## this includes cutadapt
    module load Nextflow R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0 cutadapt/3.4-GCCcore-10.3.0

    ##this includes git for interactive jobs (one-line)
    module load git/2.21.0-GCCcore-8.2.0-nodocs && module load Nextflow R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0 
    
    # pull and run latest pipeline (from project directory)
    git pull && nextflow run .

    # run pipeline with "dual" test data (default is "single")
    git pull && nextflow run . -resume --data_folder test_data/dual

    # run pipeline in preview mode with DAG output
    git pull && nextflow run . --data_folder test_data/dual -preview -with-dag dag.html

    ## run on a complete dataset to check how it copes with big sample # and sizes
    # copy files
    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/tephritid_metabarcoding/data /group/pathogens/IAWS/Personal/JackS/nextflow_tests/piperline_nextflow/test_data/full_teph
    # remove other flow cells
    rm -r /group/pathogens/IAWS/Personal/JackS/nextflow_tests/piperline_nextflow/test_data/full_teph/K3DVL
    rm -r /group/pathogens/IAWS/Personal/JackS/nextflow_tests/piperline_nextflow/test_data/full_teph/KMLK4
    # run pipeline
    git pull && nextflow run . --data_folder test_data/full_teph


    R

    > source("./jack_notes/.Rprofile")

Add above code to the "R sourcing" code block at the start of every R script (could I wrap this in an R script?)


TODO: check for all files previously checked in the `parameters_setup.R` script, which no longer checks because Nextflow breaks the relative paths

