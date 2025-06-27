# Usage

## Installation

Running Freyr requires the following software installed to your path:
- `bash` version `3.2` (or later)
- `nextflow` (the ['self-install'](https://www.nextflow.io/docs/latest/install.html#self-install) package)
- `java` version `17` (or later)
- `git`
- One of the following container managers: [`docker`](https://www.docker.com/), [`apptainer`](https://apptainer.org/), [`singularity`](https://sylabs.io/singularity/), [`podman`](https://podman.io/) or [`shifter`](https://docs.nersc.gov/development/containers/shifter/)


## Quick start

```
analysis_dir=/path/to/my/analysis/dir

# clone Freyr repository

git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $analysis_dir 

cd $analysis_dir

# set nextflow version

export NXF_VER=23.05.0-edge 

# run pipeline using your particular container manager and test data with a local executor

nextflow run . -profile test,[my_container_manager] -executor local  
```

To get a list of allowed parameters/commands while in your analysis directory, run:

    nextflow run . --help

## Input files

The pipeline has two main input files: the **samplesheet** file and the **primer parameters** file. Both must be comma-delimited (`.csv`) files, and are passed to the pipeline using `--samplesheet` and `--primer_params`, respectively.

### Samplesheet file

The samplesheet contains information on each sample to be analysed in a specific run of the pipeline. It can contain the following fields/headers/columns:

| Field | Necessity | Description | Specifications | Example value |
| --- | --- | --- | --- | --- |
| `sample` | Required | Name/ID of the sample | Must be unique. If `read_dir` is used, must be present at the start of the read file names. | `sampleA1` |
| `read_group` | Optional | Read group/flowcell/sequencing lane of the sample | If unsure, omit this column. | `KMX725` |
| `primers` | Required | Name of the PCR primer pair(s) used to generate amplicons for this sample | Can be a list separated by semi-colons (`;`) | `COIfwd-rev;marker3` |
| `read_dir` | Required\* | Directory containing the reads for the sample. | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads` |
| `fwd` | Required\* | Forward read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725_R1.fq.gz` |
| `rev` | Required\* | Reverse read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725_R2.fq.gz` |
| `single` | Required\* | Read file path | Must be a path, absolute or relative to the analysis directory. | `/path/to/reads/sampleA1_KMX725.fq.gz` |
| \<other\> | Optional | Additional, arbitrary metadata | Can be (nearly) any field name, with (nearly) any values. | `control`, `treatment1` |

*\* Only one of the three following column combinations should be used: `read_dir`, OR `fwd` and `rev`, OR `single`.*

> [!NOTE]
> Samplesheet fields and values may not contain commas, spaces or any other whitespace characters.

#### Sample

The easiest way to guarantee `sample` value uniqueness is to put the name of the sequencing flowcell (or read group) at the start, meaning that even if the same biological sample is sequenced twice, it will be treated as a different sample by the pipeline.

#### Read group

["Read group"](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) is used to describe a set of reads that all derive from the same lane of a sequencing flow cell, and is used in the pipeline for error correction. If you know this information, you can supply it using `read_group`, but it is usually safer to allow the pipeline to determine it from each sample's read headers. To do the latter, omit the `read_group` field from the samplesheet entirely. 

#### Primers

Each sample represents a collection of sequenced amplicons produced by PCR reactions. Freyr supports multiplexed sequencing designs, where multiple types of amplicons, derived from different PCR primer pairs, can be pooled into each sample. The `primers` field tells Freyr which primers were used to generate amplicons in each sample -- if multiple primers were used, supply them as a list delimited with `;`. Each sample can have any combination of primers, they don't have to be the same in every sample.

The values in `primers` must match (or be a subset) of those in the [primer parameters file](#primer-parameters-file).

#### Read files

You have two options when supplying read file information in the samplesheet: you can specify a directory where the reads for that sample can be found with `read_dir`, or you can give direct file paths per sample, with `fwd` and `rev` if reads are paired-end or with `single` if reads are single-end.

If `read_dir` is used, the `sample` value must be present at the start of the read file names. The `read_dir` directory is searched recursively, so if a read file is at `/home/my_data/raw_data/sampleX.fq`, you could use `/home/my_data`, as long as a sister directory like `/home/my_data/processed_data` didn't also contain a matching file.

#### Metadata

The samplesheet file can also contain additional metadata columns, which are ignored by the pipeline until the generation of output files. These can be useful to pass through additional information for each sample to make post-pipeline analyses simpler. For example, you could have a `location` field with the values `Melbourne`, `Sydney` and `Brisbane`.

Metadata columns can have (nearly) any name, as long as they're not one of the essential fields in the table [above](#samplesheet-file).  

### Primer parameters file

The primer parameters file tells the pipeline how to process amplicons derived from a particular primer pair. It can contain the following fields/headers/columns:

| Field | Necessity | Description | Specifications | Default value |
| --- | --- | --- | --- | --- |
| `primers` | Required | Name of primer pair | Must be unique within file | No default, must be supplied |
| `locus` | Required | Target locus/gene for primers | None (does not have to be unique) | No default, must be supplied |
| `for_primer_seq` | Required | 5'-3' nucleotide sequence of the forward primer | [IUPAC characters](https://genome.ucsc.edu/goldenPath/help/iupac.html) only, of any case | No default, must be supplied |
| `rev_primer_seq` | Required | 5'-3' nucleotide sequence of the reverse primer | [IUPAC characters](https://genome.ucsc.edu/goldenPath/help/iupac.html) only, of any case | No default, must be supplied |
| `ref_fasta` | Required | A reference database of taxonomically classified sequences | A path, absolute or relative to the analysis directory | No default, must be supplied |
| `max_primer_mismatch` | Optional | Maximum number of mismatches to allow when detecting primer sequences | Integer >= `0` | `0` |
| `read_min_length` | Optional | Minimum allowed length of reads after primer-trimming | Integer >= `0` | `20` |
| `read_max_length` | Optional | Maximum allowed length of reads after primer-trimming | Integer >= `0` or `Inf` | `Inf` |
| `read_max_ee` | Optional | Maximum expected errors allowed in primer-trimmed reads | Integer >= `0` | `1` |
| `read_trunc_length` | Optional | Length to which all primer-reads will be truncated (right-trimming) | Integer >= `1`, or `0` to disable | `0` |
| `read_trim_left` | Optional | Bases to be trimmed from the left side of all primer-trimmed reads | Integer >= `0` | `0` |
| `read_trim_right` | Optional | Bases to be trimmed from the right side of all primer-trimmed reads | Integer >= `0` | `0` |
| `asv_min_length` | Optional | Minimum allowable length of each ASV | Integer >= `0` | `0` |
| `asv_max_length` | Optional | Maximum allowable length of each ASV | Integer >= `0` or `Inf` | `Inf` |
| `concat_unmerged` | Optional | Retain unmerged read pairs by concatenating  separated by a string of 10 `N` bases | `TRUE`/`T` or `FALSE`/`F` | `FALSE` |
| `genetic_code` | Optional | Genetic code for amplicon, if a coding sequence | A value from [`Biostrings::GENETIC_CODE_TABLE`](https://rdrr.io/bioc/Biostrings/man/GENETIC_CODE.html), or `NA` | `NA` |
| `coding` | Optional | Whether the amplicon is a coding sequence or not | `TRUE`/`T` or `FALSE`/`F` | `FALSE` |
| `phmm` | Optional | A [profile Hidden Markov Model (PHMM)](https://github.com/shaunpwilkinson/aphid) of the amplicon | A path, absolute or relative to the analysis directory, or `NA` | `NA` |
| `idtaxa_db` | Optional | A trained [`IDTAXA` model](https://rdrr.io/bioc/DECIPHER/man/LearnTaxa.html) of `ref_fasta` database | A path, absolute or relative to the analysis directory, or `NA` | `NA` |
| `idtaxa_confidence` | Optional | Minimum bootstrap confidence for `IDTAXA` classification | Integer `0`-`100` | `60` |
| `run_blast` | Optional | Whether to run `BLAST` to complement `IDATAXA` classification | `TRUE`/`T` or `FALSE`/`F` | `FALSE` |
| `blast_min_identity` | Optional | Minimum nucleotide identity (%) for `BLAST` taxonomic assignment | Integer `0`-`100` | `97` |
| `blast_min_coverage` | Optional | Minimum query coverage (%) for `BLAST` taxonomic assignment | Integer `0`-`100` | `90` |
| `cluster_threshold` | Optional | ASV similarity clustering threshold (%) | Integer `0`-`100`, or `NA` | `NA` |
| `target_kingdom` | Optional | Filter output ASVs to this taxonomic kingdom | A single valid taxon name | `NA` |
| `target_phylum` | Optional | Filter output ASVs to this taxonomic phylum | A single valid taxon name | `NA` |
| `target_class` | Optional | Filter output ASVs to this taxonomic class | A single valid taxon name | `NA` |
| `target_order` | Optional | Filter output ASVs to this taxonomic order | A single valid taxon name | `NA` |
| `target_family` | Optional | Filter output ASVs to this taxonomic family | A single valid taxon name | `NA` |
| `target_genus` | Optional | Filter output ASVs to this taxonomic genus | A single valid taxon name | `NA` |
| `target_species` | Optional | Filter output ASVs to this taxonomic species | A single valid taxon name | `NA` |
| `min_sample_reads` | Optional | Minimum reads per sample to pass final filter | Integer >= `0` | `0` |
| `min_taxa_reads` | Optional | Minimum reads per ASV to pass final filter | Integer >= `0` | `0` |
| `min_taxa_ra` | Optional | Minimum relative abundance per ASV/taxon to pass final filter | Float `0`-`1`, scientific notation (eg. `1e-04`) is allowed | `0` |

If required values are not provided in the primer parameters file, they must be specified using the [override parameters](#primer-parameter-overrides).

## Pipeline parameters

All pipeline parameters (passed to `nextflow` on the command line using a double-hyphen, eg. `--parameter1`) are optional except for `--samplesheet` and `--primer_params`. Technically, `--primer_params` *is* optional, but users must instead use a specific set of `--pp_*` [override parameters](#primer-parameter-overrides), which is only recommended for advanced users.

> [!NOTE] `nextflow` considers parameters passed with no value (eg. `nextflow run . --parameter1`) to have their value set to `true`. Unset parameters without a default have a value of `null`. You can forcibly unset a parameter by passing `null`, although this will often throw an error.

### Input options

These options are either mandatory to use or frequently required. `--miseq_internal` and `--miseq_dir` should only be used by AgVic researchers. 

| Parameter | Description | Specification | Default value |
| --- | --- | --- | --- |
| `--samplesheet` | Path to the samplesheet `.csv` file | Path, absolute or relative to the analysis directory | No default, must be specified |
| `--primer_params` | Path to the primer parameters `.csv` file | Path, absolute or relative to the analysis directory | No default, typically must be specified |
| `--seq_type` | Sequencing platform of input reads | `illumina`/`nanopore` | `illumina` |
| `--paired` | Whether reads are paired-end | Boolean (`true`/`false`) | `true` |
| `--train_idtaxa` | Whether to train an `IDTAXA` model from each reference database file | Boolean (`true`/`false`) | `false` |
| `--miseq_internal` | Whether data was generated using AgVic MiSeq platform | Boolean (`true`/`false`) | `false` |
| `--miseq_dir` | Path to AgVic MiSeq demultiplexed reads | Path, absolute or relative to the analysis directory | `null` |

### Pipeline analysis options

These options generally will not need to be changed by regular users. 

| Parameter | Description | Specification | Default value |
| --- | --- | --- | --- |
| `--primer_error_rate` | Maximum error rate when detecting primers ([`cutadapt -e`](https://cutadapt.readthedocs.io/en/v4.7/reference.html#adapter-finding-options)) | Number >= `0` | `1` |
| `--primer_n_trim` | Recognise `N` bases in reads when detecting primers ([`cutadapt --match-read-wildcards`](https://cutadapt.readthedocs.io/en/v4.7/reference.html#adapter-finding-options)) | Boolean (`true`/`false`) | `false` |
| `--high_sensitivity` | Infer ASVs with `dada2` [pseudo-pooling](https://benjjneb.github.io/dada2/pseudo.html) mode | Boolean (`true`/`false`) | `true` |
| `--dada_band_size` | Set `BAND_SIZE` parameter for `dada2::dada()` | Integer | `16` |
| `--dada_homopolymer` | Set `HOMOPOLYMER_GAP_PENALTY` parameter for `dada2::dada()` | Integer < `0` | `null` |
| `--chimera_sample_frac` | Minimum fraction of samples in which a sequence must be flagged as a chimera to be classified a chimera (see [`dada2::isBimeraDenovoTable(minSampleFraction)`](https://rdrr.io/bioc/dada2/man/isBimeraDenovoTable.html)) | Number `0`-`1` | `0.9` |
| `--chunk_taxassign` | Number of sequences to be taxonomically assigned per process | Integer > `0` | `100` |
| `--accumulation_curve` | Generate accumulation curves per sample (can be slow for large datasets) | Boolean (`true`/`false`) | `true` |
| `--merge_clusters` | Produce additional outputs where each ASV cluster is merged into a single representative sequence; value determines how taxonomy of cluster is merged | `central`, `lca`, `frequency`, `rank` or `abundance` | `null` |

`--primer_n_trim`, and increasing `--primer_error_rate`, can help to retain more reads from low-quality sequencing runs, but should only be used if necessary.

The [primer parameter](#primer-parameters-file) `cluster_threshold` can be used to cluster ASVs derived from each primer into groups (ie. OTUs). This fills in the `cluster` column in some of the pipeline's output files but still retains information for each ASV. If `--merge_clusters` is set, the pipeline will produce some additional 'clustered' outputs where each cluster of (only filtered) ASVs is merged into a single representative sequence. The value of `--merge_clusters` determines how taxonomic information is merged within each cluster:

- `central` directly uses the taxonomic assignment of the central/representative ASV
- `lca` uses the lowest common ancestor of all ASVs, with some ranks unclassified if necessary
- `frequency` uses the most common taxonomic assignment, with ties broken arbitrarily but reproducibly
- `rank` uses the taxonomic assignment with the lowest rank (eg. species assignment preferred to genus assignment), with ties broken by the `frequency` method
- `abundance` uses the taxonomic assignment of the ASV with the highest read abundance across the entire dataset

> [!IMPORTANT] Cluster assignment labels (ie. `1, 2, 3...`) should not be compared between separate pipeline runs. They will also be different between the 'unfiltered' and 'filtered' outputs of Freyr. If unsure, we recommend users perform their own OTU clustering.

### Primer parameter overrides

Values in the primer parameters `.csv` file (given to `--primer_params`) can be overridden on the command line with pipeline parameters of the form `--pp_<name of primer parameter>`, eg. `--pp_max_primer_mismatch`. This allows users to easily change primer parameters between runs for exploratory purposes, without needing to modify the `--primer_params` file.

Override parameters can be used in a various ways. The simplest is to give a single value, which will be used to override that field for all rows in the `--primer_params` file, eg. `--pp_read_max_ee 2`.

To replace the values of one or more specific rows, you can 'tag' values with the primers name in the form `[<primer_name>]<parameter_value>`. For example, to set `blast_min_identity` to `95` for the primer pair `COIfwd-rev` only, you would use `--pp_blast_min_identity [COIfwd-rev]95`.

Multiple specific values can be given by separating with `;`, although the value must be wrapped in quotes, eg. `--pp_read_min_length "[primersA]30;[primersB]40"`. Primers not tagged will not be overridden.

If `--primer_params` is unset, you can still run the pipeline as long as `--pp_primers`, `--pp_locus`, `--pp_for_primer_seq`, `--pp_rev_primer_seq` and `--pp_ref_fasta` are all set appropriately. However, we strongly recommend the use of `--primer_params` to simplify pipeline execution, especially when more than one primer pair is used. To prevent complications, `--pp_primers` cannot be used unless `--primer_params` is unset.

This can all get very complicated -- when in doubt, just use a `--primer_params` file!

### Specifying resource limits

| Parameter | Description | Specification | Default value |
| --- | --- | --- | --- |
| `--max_memory` | Maximum memory available for any individual process | Of the form `<number>.<MemoryUnit>` (see [here](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#memoryunit)) | `128.GB` |
| `--max_cpus` | Maximum CPUs available for any individual process | Integer >= `1` | `16` |
| `--max_time` | Maximum time available for any individual process | Of the form `<number>.<DurationUnit>` (see [here](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#duration)) | `240.h` |

You should be careful to set `--max_memory`,`--max_cpus` and/or `--max_time` appropriately for your computational environment. For example, if you're using a laptop with 64GB RAM and 16 CPU cores, you should set `--max_memory` to less than 64GB and `--max_cpus` to less than 16. In some cases, setting these too low may cause the pipeline to fail.

### Debugging options

These parameters can help with debugging or configuring the pipeline.

| Parameter | Description | Specification | Default value |
| --- | --- | --- | --- |
| `--rdata` | Saves an `.rda` RData file in the work directory of each R-scripted process, regardless of exit status | Boolean (`true`/`false`) | `false` |
| `--subsample` | Pseudorandomly reduce number of input samples per `read_group` x `primers` combination to this | Integer >= `1` | `null` |
| `--downsample` | Pseudorandomly reduce number of reads per input sample to this | Integer >= `1` | `null` |

## Profiles

Nextflow uses [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to set collections of pipeline parameters all at once. This is useful to configure the pipeline for particular running situations (eg. cluster vs. laptop, real data vs. test data). Profiles are defined on the command line with the `-profile` flag (note the single hyphen). Multiple profiles can be used at once, separated by commas, but their ordering matters: later profiles override the settings of earlier profiles.

For example, to use both the `basc_slurm` profile (for running on BASC with the SLURM executor) and `test` profile (for running a minimal test dataset included with the pipeline), you would specify `-profile basc_slurm,test` when running the pipeline. Because `test` comes second, it overrides the max job request parameters (eg. `params.max_memory`) specified by `basc_slurm`, which is useful in this case because it will likely make job allocation through SLURM much faster.

Profiles for using [`docker`](https://www.docker.com/), [`apptainer`](https://apptainer.org/), [`singularity`](https://sylabs.io/singularity/), [`podman`](https://podman.io/) and [`shifter`](https://docs.nersc.gov/development/containers/shifter/) are all available, eg. `-profile shifter`. 

You can create and use custom profiles by [writing your own](https://www.nextflow.io/docs/latest/config.html) Nextflow `.config` file and specifying it with `-c path/to/config/file` when running `freyr`. A tutorial on how to do this will be available soon.
