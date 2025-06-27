# `freyr`: a metabarcoding analysis pipeline for agricultural biosecurity and biosurveillance

<center><img src="./assets/images/freyr.png" alt="The god Freyr, riding his boar, Gullinbursti: Murray, Alexander (1874). Manual of Mythology : Greek and Roman, Norse, and Old German, Hindoo and Egyptian Mythology. London, Asher and Co. https://commons.wikimedia.org/wiki/File:Freyr_riding_Gullinbursti.jpg" width="400"/></center>

**`freyr`** is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based metabarcoding analysis pipeline, primarily designed for use in biosecurity and biosurveillance in agriculture. It intends to allow highly reproducible, scalable, user-friendly and interpetable analyses, as well as flexibility across a wide variety of metabarcoding experiments. 

This pipeline is being developed by [Agriculture Victoria Research](https://agriculture.vic.gov.au/), as a part of the [NGDSI project](https://grdc.com.au/grdc-investments/investments/investment?code=DEE2305-004RTX), with primary development by [Jack Scanlan](https://github.com/jackscanlan) and [Alex Piper](https://github.com/alexpiper). `freyr` is the successor to [`pipeRline`](https://github.com/alexpiper/piperline).

### News

**2025-06-19:** The samplesheet and primer parameter inputs to `freyr` have recently changed -- please see the [Usage](/docs/usage.md) page for detailed information on how to run the pipeline. 

**2024-08-26:** A [step-by-step guide/tutorial](/docs/nanopore_tutorial.md) focused on analyses of Nanopore data is now available for AgVic users of the pipeline with access to the BASC HPC system.

**2024-08-09:** A [step-by-step guide/tutorial](/docs/insect_coi.md) focused on typical insect COI analyses is now available for AgVic users of the pipeline with access to the BASC HPC system. 


### Introduction

> This pipeline is currently **experimental** and being actively developed, with no stable releases as yet! If you need a stable metabarcoding pipeline, we currently recommend [`pipeRline`](https://github.com/alexpiper/piperline) or [`nfcore/ampliseq`](https://github.com/nf-core/ampliseq).

`freyr` is designed for analysing amplicon-based metabarcoding data where sequencing libraries were prepared without fragmentation. Currently short-read (eg. Illumina, MGI) paired- and single-end data are best supported, but long-read Nanopore data is supported in a very experimental state. 

This pipeline currently requires the use of `nextflow` version `23.05.0-edge` on a Linux AMD64/x86-64 system, and one of the following container managers: [`docker`](https://www.docker.com/), [`apptainer`](https://apptainer.org/), [`singularity`](https://sylabs.io/singularity/), [`podman`](https://podman.io/) or [`shifter`](https://docs.nersc.gov/development/containers/shifter/).  

### Quick start

Install a [self-install `nextflow` package](https://www.nextflow.io/docs/latest/install.html#self-install) and then run the following:

```
analysis_dir=/path/to/my/analysis/dir

# clone Freyr repository
git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $analysis_dir 

cd $analysis_dir

# set nextflow version
export NXF_VER=23.05.0-edge 

# run pipeline
nextflow run . \
    --samplesheet my_samplesheet.csv \
    --primer_params my_primer_params.csv \
    -profile shifter    
```

To get a list of allowed parameters/commands while in your analysis directory, run:

    nextflow run . --help

For detailed information on how to run the pipeline, see [Usage](/docs/usage.md). 

### Development Roadmap

The following are earmarked for inclusion into the pipeline in the near future. If you would like to request a specific feature, please [open an issue](https://github.com/AVR-biosecurity-bioinformatics/freyr/issues).
- Support for additional taxonomic classifiers other than `IDTAXA` and `BLAST`
- Support for newer versions of `nextflow`
- Increased barcode flexibility
- Improved reporting outputs for biosecurity and biosurveillance applications
- Improved support for Nanopore data

