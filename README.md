# PipeRline

`pipeRline` is a `targets`-based metabarcoding pipeline currently implemented in R. 

This pipeline is currently being re-written in the **Nextflow** language as a _containerised_ and _reproducible_ workflow. It is also partially inspired by the [nfcore/ampliseq](https://github.com/nf-core/ampliseq) pipeline!


### Usage

**2024-02-27:** The `R:4.2.0` environment (and `blast+`) is now available in a (linux/amd64) Docker image available at [Docker Hub](https://hub.docker.com/repository/docker/jackscanlan/piperline/general). Scripts are provided to run the pipeline on the BASC HPC, which uses [SLURM](https://slurm.schedmd.com/) and [Shifter](https://github.com/NERSC/shifter), in the `running_scripts` directory. 

To run the pipeline in a SLURM-based HPC using `shifter`, modify the `RUN_LOC` and `RUN_DIR` variables in `basc_shifter.slurm`, and also modify the `module` commands to be specific for your system. 







## old README text

The best place to start is going through the [General introduction to the pipeRline workflow](https://alexpiper.github.io/piperline/vignettes/general.html)

Project specific workflows:

* [Insect COI](https://alexpiper.github.io/piperline/vignettes/insect_coi.html)
 
* [Tephritid surveillance (COI + EIF3L)](https://alexpiper.github.io/piperline/vignettes/tephritid.html)

* [Marine surveillance (COI + 18S)](https://alexpiper.github.io/piperline/vignettes/marine_surveillance.html)

* [Bee health (16S + ITS)](https://alexpiper.github.io/piperline/vignettes/fungal_bacterial.html)
