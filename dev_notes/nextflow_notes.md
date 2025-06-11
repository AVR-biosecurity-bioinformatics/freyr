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
    module load git/2.21.0-GCCcore-8.2.0-nodocs && module load Nextflow/24.04.1 R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0 
    
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

    # development of branch
    git pull origin direct_package_loading && nextflow run . -resume --data_folder test_data/dual

    R

    > source("./jack_notes/.Rprofile")

Add above code to the "R sourcing" code block at the start of every R script (could I wrap this in an R script?)


TODO: check for all files previously checked in the `parameters_setup.R` script, which no longer checks because Nextflow breaks the relative paths

### setting up `/home/js7t/personal/dev/piperline`, home of development for local code testing pre-commit 

    # make dir
    mkdir -p /home/js7t/personal/dev/piperline && cd /home/js7t/personal/dev/piperline

    # clone repo
    git clone https://github.com/jackscanlan/piperline.git .

    # copy reference databases
    cp -r /home/js7t/personal/nextflow_tests/piperline_nextflow/reference/* /home/js7t/personal/dev/piperline/reference

    # copy full_teph dataset
    cp -r /home/js7t/personal/nextflow_tests/piperline_nextflow/test_data/full_teph /home/js7t/personal/dev/piperline/test_data

Running code:

    # go to folder 
    cd /home/js7t/personal/dev/piperline
    
    # interactive session
    sinteractive -c16
    
    # load modules
    module load git/2.21.0-GCCcore-8.2.0-nodocs && module load Nextflow R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0

    # remove old outputs
    rm -rf work/* output/modules/*

    # run tests_data/dual
    nextflow run . -resume --data_folder test_data/dual


### installing Nextflow for personal use in BASC

    cd ~
    module load Java/17.0.6                     # need Java v11+ to use Nextflow
    curl -s https://get.nextflow.io | bash      # install nextflow in currect dir
    chmod 777 nextflow                          # make executable
    # mkdir -p /home/js7t/bin                   # make local bin
    mv nextflow /home/js7t/bin                  # move into executable path

Running Nextflow with local version:

    module load Java/17.0.6 R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0

    nextflow run . -resume --data_folder test_data/dual

    nextflow run . -resume --samplesheet test_data/dual/samplesheet_read_dir.csv


Run with shifter:

    module load Java/17.0.6 shifter/22.02.1

    NXF_VER=23.05.0-edge \
        nextflow run . \
        --samplesheet test_data/dual/samplesheet_read_dir.csv \
        -profile shifter

Run with shifter and slurm in BASC

    module load Java/17.0.6 shifter/22.02.1

    export NXF_EXECUTOR=slurm

    NXF_VER=23.05.0-edge \
        nextflow run . \
        --samplesheet test_data/dual/samplesheet_read_dir.csv \
        -profile basc_slurm

The above code runs the pipeline with the BASC slurm system, such that each process queues a new job. 
- TODO: add better process options (cpu, time, memory) to each module, rather than blanket ones across the whole pipeline
    - use dynamic directives for this: https://www.nextflow.io/docs/latest/process.html#dynamic-computing-resources

Run with shifter and slurm in BASC, with test data

    module load Java/17.0.6 shifter/22.02.1

    NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,test,debug

Run test data with minimal samples:

    module load Java/17

    NXF_VER=23.05.0-edge nextflow run . -profile basc_slurm,test,debug --subsample 1



Test with Horsham data (one flowcell, three primer pairs)

    nextflow run . \
        -profile basc_slurm,debug \
        --samplesheet /group/pathogens/IAWS/Projects/NGDSI/Horsham/samplesheets/freyr_samplesheet_LM2TV.csv \
        --loci_params /group/pathogens/IAWS/Projects/NGDSI/Horsham/analysis/loci_params.csv \
        --lp_idtaxa_db /group/pathogens/IAWS/Personal/JackS/databases/coi/tarth_250307/idtaxa_model.rds \
        --lp_ref_fasta /group/pathogens/IAWS/Personal/JackS/databases/coi/tarth_250307/final_database.fasta \
        -resume