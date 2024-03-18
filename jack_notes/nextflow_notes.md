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
    module load Nextflow R/4.2.0-foss-2021b pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2 GDAL/3.3.0-foss-2021a BLAST+/2.11.0-gompi-2020a Pandoc/2.5 ZeroMQ/4.3.2-GCCcore-9.3.0 cutadapt/3.4-GCCcore-10.3.0
    ## this includes cutadapt

    # pull and run latest pipeline
    git pull && nextflow run .

    R

    > source("./jack_notes/.Rprofile")

Add above code to the "R sourcing" code block at the start of every R script (could I wrap this in an R script?)


TODO: check for all files previously checked in the `parameters_setup.R` script, which no longer checks because Nextflow breaks the relative paths

