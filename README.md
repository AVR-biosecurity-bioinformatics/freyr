# `freyr`: a metabarcoding analysis pipeline for agricultural biosecurity and biosurveillance

<center><img src="./assets/images/freyr.png" alt="The god Freyr, riding his boar, Gullinbursti: Murray, Alexander (1874). Manual of Mythology : Greek and Roman, Norse, and Old German, Hindoo and Egyptian Mythology. London, Asher and Co. https://commons.wikimedia.org/wiki/File:Freyr_riding_Gullinbursti.jpg" width="400"/></center>

`freyr` is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based metabarcoding analysis pipeline, primarily designed for use in biosecurity and biosurveillance in agriculture. It is the successor to [`pipeRline`](https://github.com/alexpiper/piperline) and is also inspired by [`nfcore/ampliseq`](https://github.com/nf-core/ampliseq). `freyr` intends to be highly reproducible, scalable, user-friendly and interpetable, as well as flexible across a wide variety of metabarcoding experiments. 

This pipeline is being developed by a team at [Agriculture Victoria Research](https://agriculture.vic.gov.au/), as a part of the [National Grain Diagnostic & Surveillance Initiative (NGDSI)](https://grdc.com.au/grdc-investments/investments/investment?code=DEE2305-004RTX). 


### Usage

> This pipeline is currently **experimental** and being actively developed, with no guarantee that the code is stable! If you need a stable metabarcoding pipeline, we currently recommend [`pipeRline`](https://github.com/alexpiper/piperline).

Running `freyr` might look something like this:

    nextflow run AVR-biosecurity-bioinformatics/freyr \
        --samplesheet samplesheet.csv \
        --loci_params loci_params.csv \
        -profile test

**2024-02-27:** The `R:4.2.0` environment (and `blast+`) is now available in a (linux/amd64) Docker image available at [Docker Hub](https://hub.docker.com/repository/docker/jackscanlan/piperline/general). Scripts are provided to run the pipeline on the BASC HPC, which uses [SLURM](https://slurm.schedmd.com/) and [Shifter](https://github.com/NERSC/shifter), in the `running_scripts` directory. 

**2024-06-12:** To run the pipeline completely using containers (which is recommended), install [`nextflow`](https://www.nextflow.io/docs/latest/install.html) (not the `all` distribution, in order to allow third-party plugins) and [`shifter`](https://shifter.readthedocs.io/en/latest/install_guides.html). Use the following commands to run a test dataset:

    ## make sure java, nextflow and shifter are all available in your path
    
    # run the pipeline using nextflow v23.04.5
    NXF_VER=23.04.5 \
        nextflow run AVR-biosecurity-bioinformatics/freyr \
        --samplesheet test_data/dual/samplesheet_read_dir.csv \
        -profile shifter

The pipeline may also work with `-profile` set to `apptainer`, `conda`, `docker`, `podman` or `singularity`--when using the respective container platform--but these have currently not been tested internally. 


#### Important notes:

- This pipeline currently only works with native Shifter support (ie. with `-profile shifter` in the Nextflow run command) if Nextflow is version `23.04.5` (or possibly older). This is due to a bug in how Nextflow (at least versions `23.10.0` to `24.04.2`) sets up the process environment in `.command.run`
- The pipeline has not been tested with Docker, Singularity, Apptainer or Podman--only Shifter. If you run the pipeline using one of these platforms, please let us know if it works or not!
- Charliecloud has been tested and will not work with this pipeline. 
    - It seems like Nextflow (as of version `24.04.2`) fails to handle using Charliecloud with pipelines that use more than one container, as they cannot all get pulled at the same time. Even for pipelines that use only a single container, moving the container contents to each process work directory slows the execution down considerably and increases the project footprint. 
- When running the pipeline with containers, you must be using a Linux system with AMD64 architecture (such as AgVic's BASC). In the future, we will try to support other architectures by using multi-platform containers. 


### Inputs

The pipeline has two main inputs: the **samplesheet**, and the **loci parameters**.  

The samplesheet tells the pipeline what samples are being run, as well as (for each sample): where the sequencing read files are, what primers were used, what flowcell/experiment they were sequenced in, and additional metadata. The samplesheet should be provided to the pipeline as a `.csv` file using the `--samplesheet` flag, where each row is a different sample.

The loci parameters tell the pipeline how to analyse the samples, on a per-locus basis (for multiplexed experiments where multiple loci were pooled per sample). The loci parameters should be provided to the pipeline as a `.csv` file using the `--loci_params` flag, where each row is a different locus/PCR primer pair. 

Both samplesheet and loci parameters `.csv` files are checked by the pipeline before the run starts, to make sure all the values provided are valid and what the pipeline will expect. 


### General parameters

#### Resource limits

If your computational environment has hard limits on the resources it can devote to the pipeline (eg. you're running on a personal computer with a relatively small amount of CPU and memory), you should be careful to set `params.max_memory`,`params.max_cpus` and/or `params.max_time`. This will make sure the pipeline as a whole (for local execution), or any particular process (for cluster/SLURM execution), stays within these limits.  

By default these are set to:
- `params.max_memory = 128.GB`
- `params.max_cpus = 16`
- `params.max_time = 240.h`

### Profiles

Nextflow uses [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to set collections of pipeline parameters all at once. This is useful to configure the pipeline for particular running situations (eg. cluster vs. laptop, real data vs. test data). Profiles are defined on the command line with the `-profile` flag. Multiple profiles can be used at once, separated by commas, but their ordering matters: later profiles override the settings of earlier profiles. 

For example, to use both the `basc_slurm` profile (for running on BASC with the SLURM executor) and `test` profile (for running a minimal test dataset included with the pipeline), you would specify `-profile basc_slurm,test` when running the pipeline. Because `test` comes second, it overrides the max job request parameters (eg. `params.max_memory`) specified by `basc_slurm`, which is useful in this case because it will likely make job allocation through SLURM much faster.

You can create and use custom profiles by [writing your own](https://www.nextflow.io/docs/latest/config.html) Nextflow `.config` file and specifying it with `-c path/to/config/file` when running `freyr`. A tutorial on how to do this will be available soon.


## old README text for `pipeRline`

The best place to start is going through the [General introduction to the pipeRline workflow](https://alexpiper.github.io/piperline/vignettes/general.html)

Project specific workflows:

* [Insect COI](https://alexpiper.github.io/piperline/vignettes/insect_coi.html)
 
* [Tephritid surveillance (COI + EIF3L)](https://alexpiper.github.io/piperline/vignettes/tephritid.html)

* [Marine surveillance (COI + 18S)](https://alexpiper.github.io/piperline/vignettes/marine_surveillance.html)

* [Bee health (16S + ITS)](https://alexpiper.github.io/piperline/vignettes/fungal_bacterial.html)
