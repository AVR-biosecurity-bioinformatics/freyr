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




## old README text

The best place to start is going through the [General introduction to the pipeRline workflow](https://alexpiper.github.io/piperline/vignettes/general.html)

Project specific workflows:

* [Insect COI](https://alexpiper.github.io/piperline/vignettes/insect_coi.html)
 
* [Tephritid surveillance (COI + EIF3L)](https://alexpiper.github.io/piperline/vignettes/tephritid.html)

* [Marine surveillance (COI + 18S)](https://alexpiper.github.io/piperline/vignettes/marine_surveillance.html)

* [Bee health (16S + ITS)](https://alexpiper.github.io/piperline/vignettes/fungal_bacterial.html)
