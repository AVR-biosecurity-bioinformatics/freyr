# Insect COI workflow
#### by Jack Scanlan, 2024-08-02

This tutorial will help you run `freyr` on the Agriculture Victoria BASC HPC system, for experiments utilising the 'freshwater' COI primers `fwhF2-fwhR2n` on mixed arthropod samples (focusing on insects).

## Prior to analysis

The following are some assumptions that you need to make sure are met before you start this analysis:
1. Samples have been sequenced with Illumina paired-end sequencing (eg. MiSeq or NovaSeq)
2. You have used the `fwhF2-fwhR2n` primers to amplify the COI gene, with Illumina adapters added to the end of each amplicon, ie. amplicons have **not** been fragmented before adapters were added
3. Your data is available, per sample, as gzipped FastQ files (with `.fastq.gz` or `.fq.gz` extensions), with forward and reverse reads in separate files, ie. **not** interleaved

### Install `nextflow` in your home directory

`freyr` currently requires a specific version of `nextflow` that is not available as a BASC module. It also uses a third-party validation plugin that doesn't work with standalone (ie. module-based) distributions of `nextflow`. Luckily, it is very easy to install `nextflow` for a specific user on BASC: 

```
## this can be done in a login node

# change to your home directory
cd ~

# load Java
module load Java/17.0.6

# install nextflow in current directory
curl -s https://get.nextflow.io | bash

# make nextflow executable
chmod 777 nextflow

# make a home bin directory if you don't already have one
mkdir -p ~/bin

# move nextflow into bin to make it executable from any path
mv nextflow ~/bin
```

> NOTE: This only needs to be done once before any particular user uses `nextflow` for the first time -- you don't need to repeat this step for subsequent runs of this pipeline, or any other `nextflow` pipeline.

### Clone the `freyr` GitHub repository

Each run of `freyr` is conducted in a dedicated working directory created by cloning the GitHub repository:

```
## this would usually be done in a login node
## if doing in a computational node, you need to load a `git` module first

# define working directory (example only; change to your own)
working_dir=/group/pathogens/IAWS/Personal/JackS/dev/freyr

# clone repository
git clone https://github.com/alexpiper/piperline.git $working_dir

# change to working directory
cd $working_dir
```

### Make your samplesheet

The samplesheet is a `.csv` file that contains all the information 



