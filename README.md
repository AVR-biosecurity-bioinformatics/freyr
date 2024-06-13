# PipeRline

`pipeRline` is a metabarcoding analysis pipeline largely written in R. 

This pipeline is currently being re-written in the **Nextflow** language so it will be completely _containerised_ and _reproducible_. It is also partially inspired by the [nfcore/ampliseq](https://github.com/nf-core/ampliseq) pipeline!


### Usage

**2024-02-27:** The `R:4.2.0` environment (and `blast+`) is now available in a (linux/amd64) Docker image available at [Docker Hub](https://hub.docker.com/repository/docker/jackscanlan/piperline/general). Scripts are provided to run the pipeline on the BASC HPC, which uses [SLURM](https://slurm.schedmd.com/) and [Shifter](https://github.com/NERSC/shifter), in the `running_scripts` directory. 


**2024-06-12:** To run the pipeline completely using containers (which is recommended), install [`nextflow`](https://www.nextflow.io/docs/latest/install.html) (not the `all` distribution, in order to allow third-party plugins), [`git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [`shifter`](https://shifter.readthedocs.io/en/latest/install_guides.html). Use the following commands to run a test dataset:

    # make sure git, nextflow and shifter are all available in your path
    
    # clone this GitHub project to a new directory
    git clone https://github.com/jackscanlan/piperline.git test_data && \
        cd test_data  

    # run the pipeline using nextflow v23.04.5
    NXF_VER=23.04.5 \
        nextflow run . \
        --samplesheet test_data/dual/samplesheet_read_dir.csv \
        -profile shifter

The pipeline may also work with `-profile` set to `apptainer`, `conda`, `docker`, `podman` or `singularity`--when using the respective container platform--but they have not been tested with the current code. 


#### Important notes:

- This pipeline currently only works with native Shifter support (ie. with `-profile shifter` in the Nextflow run command) if Nextflow is version `23.04.5` (or possibly older). This is due to a bug in how Nextflow (at least versions `23.10.0` to `24.04.2`) sets up the process environment in `.command.run`
- The pipeline has not been tested with Docker, Singularity, Apptainer or Podman--only Shifter. If you run the pipeline using one of these platforms, please let us know if it works or not!
- Charliecloud has been tested and will not work with this pipeline. 
    - It seems like Nextflow (as of version `24.04.2`) fails to handle using Charliecloud with pipelines that use more than one container, as they cannot all get pulled at the same time. Even for pipelines that use only a single container, moving the container contents to each process work directory slows the execution down considerably and increases the project footprint. 
- When running the pipeline with containers, you must be using a Linux system with AMD64 architecture (such as AgVic's BASC). In the future, we will try to support other architectures by using multi-platform containers. 


### Inputs

The pipeline has two main inputs: the **samplesheet**, and the **loci parameters**.  

The samplesheet tells the pipeline what samples are being run, as well as, for each sample: where the sequencing read files are, what primers were used, what flowcell/experiment they were sequenced in, and additional metadata. The samplesheet should be provided to the pipeline as a `.csv` file using the `--samplesheet` flag, where each row is a different sample.

The loci parameters tell the pipeline how to analyse the samples, on a per-locus basis (for multiplexed experiments where multiple loci were pooled per sample). The loci parameters should be provided to the pipeline as a `.csv` file using the `--loci_params` flag, where each row is a different locus/PCR primer pair. 

Both samplesheet and loci parameters `.csv` files are checked by the pipeline before the run starts, to make sure all the values provided are valid and what the pipeline will expect. 


### General parameters

#### Resource limits

If your computational environment has hard limits on the resources it can devote to the pipeline (eg. you're running on a personal computer with a set number of CPU cores and memory available), you should manually set `params.max_memory`,`params.max_cpus` and/or `params.max_time`. This will make sure the pipeline as a whole (for local execution), or any particular process (for cluster/SLURM execution), stays within these limits.  

By default these are set to:
- `params.max_memory = 128.GB`
- `params.max_cpus = 16`
- `params.max_time = 240.h`


## old README text

The best place to start is going through the [General introduction to the pipeRline workflow](https://alexpiper.github.io/piperline/vignettes/general.html)

Project specific workflows:

* [Insect COI](https://alexpiper.github.io/piperline/vignettes/insect_coi.html)
 
* [Tephritid surveillance (COI + EIF3L)](https://alexpiper.github.io/piperline/vignettes/tephritid.html)

* [Marine surveillance (COI + 18S)](https://alexpiper.github.io/piperline/vignettes/marine_surveillance.html)

* [Bee health (16S + ITS)](https://alexpiper.github.io/piperline/vignettes/fungal_bacterial.html)
