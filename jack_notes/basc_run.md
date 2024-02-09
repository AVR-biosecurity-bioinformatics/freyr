# Working notes on how to run the R pipeline in BASC

This will run in the `dros_test` directory in `/group/pathogens/IAWS/Personal/JackS/piperline_tests`. 

    # Change into the main directory you wish to make the project in
    cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

    # Clone the repository from my version
    git clone https://github.com/jackscanlan/piperline.git dros_test
    cd dros_test

Copy example demultiplexed .fastq and samplesheet from another metabarcoding project:

    cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 ./data
    cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv ./data

Load interactive job, modules and R session:

    # Create new interactive SLURM session
    sinteractive --ntasks=1 --cpus-per-task=1 --mem-per-cpu=10GB --time=72:00:00

    # load modules
    module load R/4.2.0-foss-2021b
    module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
    module load GDAL/3.3.0-foss-2021a
    module load BLAST+/2.11.0-gompi-2020a
    module load Pandoc/2.5
    module load ZeroMQ/4.3.2-GCCcore-9.3.0

    # make sure default R library exists
    mkdir -p ~/R_libs/default

    # create .Rprofile if it doesn't exist already
    echo '.libPaths("~/R_libs/default")' > ~/.Rprofile
    echo 'options(renv.config.pak.enabled = TRUE)' >> ~/.Rprofile
    # do above so renv uses pak to install packages (need to load pak first)

    # Load R
    R

    # load user source file
    > source("~/.Rprofile")

I'm avoiding using the packages `crew`, `mirai` and `nanonext` as the latter has a dependency that is missing on BASC. As such, I'm using the `_targets_nocrew.R` targets file and the `renv_nocrew.lock` file to install packages. This file now also doesn't install `taxreturn`.

Install `renv` package (version `1.0.3`):

    > install.packages("renv", version = "1.0.3", repos = "http://cran.rstudio.com/")

Install `pak` package:
    > install.packages("pak", repos = "http://cran.rstudio.com/")

Load packages:

    > library(devtools) # trying this to install taxreturn manually
    > library(renv)
    > library(pak)

Install `taxreturn` manually using `devtools`:

    > install_github("alexpiper/taxreturn@e9dc03a", force = T)

Load taxreturn independently from personal library:

    > library(taxreturn)

Restore packages from `renv_nocrew.lock` file using `renv`:

    > renv::restore(lockfile = "./renv_nocrew.lock")
    # type 'y' to activate project, then "y" again to proceed installation; packages will install from download or cache (takes a while)

`pak` works really well and installs packages in parallel--much faster!

Create new `_targets_packages_nocrew.R` file that doesn't contain `crew`, `mirai` or `nanonext` and source it, installing the packages:

    > source("_targets_packages_nocrew.R")
    
    # Source ancillary functions
    source("R/functions.R")
    source("R/themes.R")

Reference database from Zenodo:

    # Download database from zenodo
    > download_zenodo(
        doi = "10.5281/zenodo.7655352",
        path = "reference"
      )
