# Usage

## Installation

Running Freyr requires the following software installed to your path:
- `bash` version `3.2` (or later)
- `nextflow` (the ['self-install'](https://www.nextflow.io/docs/latest/install.html#self-install) package)
- `java` version `17` (or later)
- `git`
- one of the following container managers: [`docker`](https://www.docker.com/), [`apptainer`](https://apptainer.org/), [`singularity`](https://sylabs.io/singularity/), [`podman`](https://podman.io/) or [`shifter`](https://docs.nersc.gov/development/containers/shifter/)



## Quick start

```
analysis_dir=/path/to/my/analysis/dir

# clone Freyr repository
git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $analysis_dir 

cd $analysis_dir

# set nextflow version
export NXF_VER=23.05.0-edge 

# run pipeline using your particular container manager
nextflow run . \
    --samplesheet my_samplesheet.csv \
    --primer_params my_primer_params.csv \
    -profile [my_container_manager]    
```

To get a list of allowed parameters/commands while in your analysis directory, run:

    nextflow run . --help


## Input files

The pipeline has two main input files: the **samplesheet** file and the **primer parameters** file. Both must be comma-delimited (`.csv`) files, and are passed to the pipeline using `--samplesheet` and `--primer_params`, respectively. 

### Samplesheet file

The samplesheet can contain the following fields/headers/columns:

| Field | Necessity | Description | Specifications | Example value |
| --- | --- | --- | --- | --- |
| `sample` | Required | Name/ID of the sample | Must be unique. If `read_dir` is used, must be present at the start of the read file names. | `sampleA1` |
| `read_group` | Optional | Read group/flowcell/sequencing lane of the sample | If unsure, omit this column. | `KMX725` |
| `primers` | Required | Name of the PCR primer pair(s) used to generate amplicons for this sample | Can be a list separated by semi-colons (`;`) | `COIfwd-rev;marker3` |
| `read_dir` | Required\* | Directory containing the reads for the sample. | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads` |
| `fwd` | Required\* | Forward read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725_R1.fq.gz` |
| `rev` | Required\* | Reverse read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725_R2.fq.gz` |
| `single` | Required\* | Read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725.fq.gz` |
| [other] | Optional | Additional metadata | Can be (nearly) any field name, with (nearly) any values. | `control`, `treatment1` |

*\* Only one of the three following column combinations should be used: `read_dir`, OR `fwd` and `rev`, OR `single`.*

> [!NOTE]
> Samplesheet fields and values may not contain commas, spaces or any other whitespace characters.

#### Sample

The easiest way to guarantee `sample` value uniqueness is to put the name of the sequencing flowcell (or read group) at the start, meaning that even if the same biological sample is sequenced twice, it will be treated as a different sample by the pipeline.

#### Read group

["Read group"](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) is used to describe a set of reads that all derive from the same lane of a sequencing flow cell, and is used in the pipeline for error correction. If you know this information, you can supply it using `read_group`, but it is usually safer to allow the pipeline to determine it from each sample's read headers. To do the latter, omit the `read_group` field from the samplesheet entirely. 

#### Primers

Each sample represents a collection of sequenced amplicons produced by PCR reactions. Freyr supports multiplexed sequencing designs, where multiple types of amplicons, derived from different PCR primer pairs, can be pooled into each sample. The `primers` field tells Freyr which primers were used to generate amplicons in each sample -- if multiple primers were used, supply them as a list delimited with `;`. Each sample can have any combination of primers, they don't have to be the same in every sample.

The values in `primers` must match (or be a subset) of those in the primer parameters file supplied with `--primer_params` (see below).

#### Read files

You have two options when supplying read file information in the samplesheet: you can specify a directory where the reads for that sample can be found with `read_dir`, or you can give direct file paths per sample, with `fwd` and `rev` if reads are paired-end or with `single` if reads are single-end.

If `read_dir` is used, the `sample` value must be present at the start of the read file names. The `read_dir` directory is searched recursively, so if a read file is at `/home/my_data/raw_data/sampleX.fq`, you could use `/home/my_data`, as long as a sister directory like `/home/my_data/processed_data` didn't also contain a matching file.

#### Metadata

The samplesheet file can also contain additional metadata columns, which are ignored by the pipeline until the generation of output files. These can be useful to pass through additional information for each sample to make post-pipeline analyses simpler. For example, you could have a `location` field with the values `Melbourne`, `Sydney` and `Brisbane`.

Metadata columns can have (nearly) any name, as long as they're not one of the essential fields in the table above.  

### Primer parameters file

The loci parameters tell the pipeline how to analyse the samples, on a per-locus basis (for multiplexed experiments where multiple loci were pooled per sample). The loci parameters should be provided to the pipeline as a `.csv` file using the `--loci_params` flag, where each row is a different locus/PCR primer pair. 

Both samplesheet and loci parameters `.csv` files are checked by the pipeline before the run starts, to make sure all the values provided are valid and what the pipeline will expect. 


## Pipeline parameters



### Resource limits

If your computational environment has hard limits on the resources it can devote to the pipeline (eg. you're running on a personal computer with a relatively small amount of CPU and memory), you should be careful to set `params.max_memory`,`params.max_cpus` and/or `params.max_time`. This will make sure the pipeline as a whole (for local execution), or any particular process (for cluster/SLURM execution), stays within these limits.  

By default these are set to:
- `params.max_memory = 128.GB`
- `params.max_cpus = 16`
- `params.max_time = 240.h`

## Profiles

Nextflow uses [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to set collections of pipeline parameters all at once. This is useful to configure the pipeline for particular running situations (eg. cluster vs. laptop, real data vs. test data). Profiles are defined on the command line with the `-profile` flag. Multiple profiles can be used at once, separated by commas, but their ordering matters: later profiles override the settings of earlier profiles. 

For example, to use both the `basc_slurm` profile (for running on BASC with the SLURM executor) and `test` profile (for running a minimal test dataset included with the pipeline), you would specify `-profile basc_slurm,test` when running the pipeline. Because `test` comes second, it overrides the max job request parameters (eg. `params.max_memory`) specified by `basc_slurm`, which is useful in this case because it will likely make job allocation through SLURM much faster.

You can create and use custom profiles by [writing your own](https://www.nextflow.io/docs/latest/config.html) Nextflow `.config` file and specifying it with `-c path/to/config/file` when running `freyr`. A tutorial on how to do this will be available soon.
