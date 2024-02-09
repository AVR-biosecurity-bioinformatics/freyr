# Working notes on how to run the R pipeline in BASC

This will run in the `dros_test` directory in `/group/pathogens/IAWS/Personal/JackS/piperline_tests`. 

    # Change into the main directory you wish to make the project in
    $ cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

    # Clone the repository from my version
    $ git clone https://github.com/jackscanlan/piperline.git dros_test
    $ cd dros_test

Copy example demultiplexed .fastq and samplesheet from another metabarcoding project:

    $ cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 .
    $ cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv .

Load interactive job, modules and R session:

    # Create new interactive SLURM session
    $ sinteractive --ntasks=1 --cpus-per-task=10 --mem-per-cpu=10GB --time=72:00:00

    # load modules
    $ module load R/4.2.0-foss-2021b
    $ module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
    $ module load GDAL/3.3.0-foss-2021a
    $ module load BLAST+/2.11.0-gompi-2020a
    $ module load Pandoc/2.5
    $ module load ZeroMQ/4.3.2-GCCcore-9.3.0
    $ module load CMake/3.26.3-GCCcore-12.3.0 # trying this to get 'nanonext' package to install correctly through renv
    # 'CMake/3.10.3-intel-2019a' didn't work
    # ' CMake 3.13 or higher is required.  You are running version 3.10.3 '

    # Load R
    $ R

I'm avoiding using the packages `crew` and `nanonext` as the latter has a dependency that is missing on BASC. As such, I'm using the `_targets_nocrew.R` targets file and the `renv_nocrew.lock` file to install packages.

Install `renv` package (version `1.0.3`) in my `piperline_nocrew` library in home directory:

    > install.packages("renv", version = "1.0.3", repos = "http://cran.rstudio.com/", lib = "~/R_libs/piperline_nocrew")

Load `renv` R package:

    > library(renv, lib.loc = "~/R_libs/4.2.0")

Restore packages from `renv.lock` file using `renv`:

    > renv::restore()
    # type 'y' to proceed; packages will install from download or cache (takes a while)

