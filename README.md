# `freyr`: a metabarcoding analysis pipeline for agricultural biosecurity and biosurveillance

<center><img src="./assets/images/freyr.png" alt="The god Freyr, riding his boar, Gullinbursti: Murray, Alexander (1874). Manual of Mythology : Greek and Roman, Norse, and Old German, Hindoo and Egyptian Mythology. London, Asher and Co. https://commons.wikimedia.org/wiki/File:Freyr_riding_Gullinbursti.jpg" width="400"/></center>

`freyr` is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based metabarcoding analysis pipeline, primarily designed for use in biosecurity and biosurveillance in agriculture. It is the successor to [`pipeRline`](https://github.com/alexpiper/piperline) and is also inspired by [`nfcore/ampliseq`](https://github.com/nf-core/ampliseq). `freyr` intends to allow highly reproducible, scalable, user-friendly and interpetable analyses, as well as flexibility across a wide variety of metabarcoding experiments. 

This pipeline is being developed by a team at [Agriculture Victoria Research](https://agriculture.vic.gov.au/), as a part of the [National Grains Diagnostic & Surveillance Initiative (NGDSI)](https://grdc.com.au/grdc-investments/investments/investment?code=DEE2305-004RTX). 


### Usage

> This pipeline is currently **experimental** and being actively developed, with no guarantee that the code is stable! If you need a stable metabarcoding pipeline, we currently recommend [`pipeRline`](https://github.com/alexpiper/piperline).

Running `freyr` might look something like this, if `nextflow` and `java` are in your path, and you use the container platform software `shifter`:

```
# clone repository into analysis directory
git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $analysis_dir \
    && cd $analysis_dir

# run pipeline
NXF_VER=23.05.0-edge \
    nextflow run . \
    --samplesheet samplesheet.csv \
    --loci_params loci_params.csv \
    -profile shifter    
```

The pipeline may also work with `-profile` set to `apptainer`, `docker`, `podman` or `singularity`--when using the respective container platform--but these have currently not been tested internally. 

To get a list of allowed parameters/command while in your analysis directory:

    nextflow run . --help

**2024-08-26:** A [step-by-step guide/tutorial](/docs/nanopore_tutorial.md) focused on analyses of Nanopore data is now available for AgVic users of the pipeline with access to the BASC HPC system.

**2024-08-09:** A [step-by-step guide/tutorial](/docs/insect_coi.md) focused on typical insect COI analyses is now available for AgVic users of the pipeline with access to the BASC HPC system. 

#### Important notes:

- `freyr` currently only works on data where sequencing adapters have been ligated onto the end of each amplicon (ie. fragmentation-based library preps are not supported)
- Short-read (eg. Illumina) paired-end data is currently best supported, but there is (very) experimental support for Nanopore data
- This pipeline currently only works with native Shifter support (ie. with `-profile shifter` in the Nextflow run command) if Nextflow is version `23.04.5` (or possibly older). This is due to a bug in how Nextflow (at least versions `23.10.0` to `24.04.2`) sets up the process environment in `.command.run`
- The pipeline has not been tested with Docker, Singularity, Apptainer or Podman--only Shifter. If you attempt to run the pipeline using one of these platforms, please let us know if it works or not!
- When running the pipeline with containers, you must be using a Linux system with AMD64/x86-64 architecture (such as AgVic's BASC). In the future, we aim to support other architectures by using multi-platform containers. 

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
