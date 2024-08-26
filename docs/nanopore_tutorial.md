# Nanopore tutorial [INCOMPLETE--DO NOT USE]
#### by Jack Scanlan, 2024-08 

This tutorial will help you run `freyr` on the Agriculture Victoria BASC HPC system, for experiments on mixed arthropod samples (focusing on insects).

> NOTE: Nanopore functionality in `freyr` is under active development and may not produce reliable results. In particular, ASV denoising and chimera detection are unlikely to be operating properly at the moment, so treat results with care. 

## Prior to analysis

The following are some assumptions that you need to make sure are met before you start this Nanopore-specific analysis:
1. Samples have been sequenced with Oxford Nanopore Technologies (ONT; 'Nanopore') sequencing (ie. from MinION, PromethION etc.).
2. Adapters sequences have been added/ligated to the end of each full-length amplicon, ie. amplicons have **not** been fragmented before adapters were added. 
    - This means libraries generated with 'Rapid Sequencing' kits (eg. SQK-RAD114) are **not** supported; please use 'Ligation Sequencing' kits (eg. SQK-LSK114) instead.
3. Your data is available, per sample, as (optionally gzipped) individual FastQ files (with `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` extensions).

--- 

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

Each run of `freyr` is conducted in a dedicated (preferably new and empty) working directory created by cloning the GitHub repository:

```
## this would usually be done in a login node
## if doing in a computational node, you need to load a `git` module first

# define working directory (example only; change to your own)
working_dir=/group/pathogens/IAWS/Personal/JackS/dev/freyr

# clone repository
git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $working_dir

# change to working directory
cd $working_dir
```

### Make your samplesheet

The samplesheet is a `.csv` file that conveys information to the pipeline about each sample. One row should be used per sample. A blank samplesheet can be found at `./inputs/samplesheet_blank.csv`. There are 8-9 required samplesheet columns (or 'fields'), listed in the table below. The rest of the fields are optional and can be used to specify additional metadata if you have it. 

| field | meaning | requirements | change to: |
| --- | --- | --- | --- |
| `sample_id` | Unique sample ID | Currently must be `fcid` and `sample_name` separated by an underscore. If `read_dir` is used over `single`, it must be present at the start of the name of each read file. | Depends on data |
| `sample_name` | Sample name | Has to be unique within a flowcell. See `sample_id`. | Depends on data |
| `fcid` | Flowcell ID | One of the two core groupings within the pipeline, along with `pcr_primers`. See `sample_id`. | Depends on data  |
| `target_gene` | Name of the target gene(s) amplified | Can specify multiple by separating with a semi-colon. | Name of your target gene (eg. `COI`)  |
| `pcr_primers` | Name of the PCR primer pair(s) | Can specify multiple by separating with a semi-colon. | Name of your primer pair(s) (eg.`fwhF2-fwhR2n`) |
| `for_primer_seq` | 5'-3' sequence of the forward primer(s) | Can specify multiple by separating with a semi-colon. | Your sequence (eg. `GGDACWGGWTGAACWGTWTAYCCHCC`)  |
| `rev_primer_seq` | 5'-3' sequence of the reverse primer(s) | Can specify multiple by separating with a semi-colon. | Your sequence (eg. `GTRATWGCHCCDGCTARWACWGG`)  |
| `read_dir` | Directory containing sequencing reads | If this is specified, the pipeline searches this directory (including subdirectories) for read files starting with `sample_id`. Cannot be used in conjunction with `single`; delete field if not using. | Up to you! |
| `single` | Exact path of the read file for this sample | Cannot be used in conjunction with `read_dir`; delete field if not using. | Up to you!  |

The easiest way to specify the reads for each sample is to copy your read files into the `./data` directory and then put `./data` in the `read_dir` field for every row in the samplesheet. 

> IMPORTANT NOTE: You must only use one of either `read_dir` OR `single` in the samplesheet. Delete the fields/columns you're not using, including the `fwd` and `rev` columns, which are used for paired-end sequencing data. 

> NOTE: Paths for the `read_dir` and `single` fields can be absolute or relative.

Once you have made your samplesheet `.csv` (it can have any name), upload it to the `inputs` directory.

### Make your loci parameters file

The loci parameters (or `loci_params`) file is a `.csv` that conveys information to the pipeline about how to process each type of amplicon found in each sample (ie. those derived from `pcr_primers` in the samplesheet). One row should be used per locus/primer pair. Paths can be absolute, or relative to the analysis directory. 

A 'default' file can be found at `./inputs/loci_params_default.csv` that contains defaults for 'typical' runs focusing on the COI locus, but you will need to change some of these:

| field | meaning | change to: |
| --- | --- | --- |
| `read_min_length` | Minimum length (bp) of reads to retain after primer trimming | A realistic minimum size of your target locus |
| `read_max_length` | Maximum length (bp) of reads to retain after primer trimming | A realistic maximum size of your target locus |
| `read_trunc_length` | Length (bp) to truncate all reads to | `0` (means disabled) |
| `asv_min_length` | Minimum allowed ASV length (bp) | Typically should be same as `read_min_length` |
| `asv_max_length` | Maximum allowed ASV length (bp) | Typically should be same as `read_max_length` |
| `phmm` | Path to PHMM model of your target sequence | eg. `folmer_fullength_model.rds`; if you don't have one, set to `NA` |
| `idtaxa_db` | Path to trained IDTAXA database | eg. `reference/idtaxa_bftrimmed.rds`; if you don't have one, set to `NA` |
| `ref_fasta` | Path to reference sequence database  | eg. `reference/insecta_hierarchial_bftrimmed.fa.gz` |

> NOTE: Because Nanopore reads tend to be highly variable in length, it is important to set `read_min_length`/`read_max_length` appropriately so that ultra-long, off-target reads get filtered out early in the pipeline. The presence of such reads can slow down the pipeline or even cause it to fail. 

> NOTE: If you don't have a trained IDTAXA model, you must have a [correctly formatted](https://alexpiper.github.io/taxreturn/articles/taxreturn.html) reference database `.fasta` file specified in `ref_fasta`, and you need to run the pipeline command with the additional flag `--train_idtaxa formatted`. 

Once you have made your loci parameters `.csv` (it can have any name), upload it to the `inputs` directory.

### Run pipeline

You should always run `freyr` on BASC in a computational node, not a login node. There are two main ways to do this: 
1. in an interactive shell using `sinteractive`, or 
2. with a SLURM script using `sbatch`. 

However, the commands to run the pipeline will remain the same for both options.  

#### Main pipeline command

To run `freyr`, you need to change to your analysis directory, load a Java module, then use the `nextflow run` command:

```
# change to your analysis directory
cd $working_dir

# load a Java module
module load Java/17.0.6

# run the pipeline
# make sure you replace the square bracketed file names (including the brackets) with the names of the files you made earlier
NXF_VER=23.04.5 \
    nextflow run . \
    --samplesheet ./inputs/[your_samplesheet].csv \
    --loci_params ./inputs/[your_loci_params].csv \
    -profile basc_slurm \
    --seq_type nanopore \
    --paired false

```

The above code is doing the following:
-  `NXF_VER=23.04.5` tells `nextflow` to run as version 23.04.5, which is currently required for the pipeline to work correctly on BASC
- `nextflow run .` tells `nextflow` to run the pipeline using files in the current directory
- `--samplesheet` and `--loci_params` specify where your samplesheet and loci parameters files are; these can be relative or absolute paths, here we use relative paths
- `-profile basc_slurm` tells `nextflow` to run using settings that make sense for BASC; in particular, it will spawn a new job for every process, making its computational execution very efficient, and will use `shifter` to run Docker containers for all the software tools
- `--seq_type nanopore` and `--paired false` tell the pipeline we're using Nanopore data

> NOTE: Confusingly, pipeline parameters like `--samplesheet` need to be typed on the command line using a double hyphen, while `nextflow` parameters like `-profile` need to be typed with a single hyphen. When in doubt, copy-paste from this page.

#### Advanced pipeline parameters

`freyr` has optional, advanced parameters that can be specified on the command line after `-profile`. Below are some relevant ones, but you typically won't need to alter these.

| parameter | meaning | requirements |
| --- | --- | --- |
| `--primer_error_rate` | Error rate allowed when detecting primer sequences | Sets the value for `cutadapt -e`; see the [`cutadapt` reference manual](https://cutadapt.readthedocs.io/en/v4.7/guide.html#error-tolerance) for more info |
| `--primer_n_trim` | Recognise N bases in reads when detecting primers | Can be `false` (default) or `true`; only use this if your sequencing data quality is very low |
| `--subsample` | Reduce the number of input samples per `pcr_primers` x `fcid` combination to this number | Should be a number >=1; this is mainly used for development and debugging and should only be used if you're trying to quickly work out if you have the right parameters for a particular dataset |
| `--downsample_reads` | Reduce the number of input reads per sample to at most this number | Should be a number >=1; this is mainly used for development and debugging and should only be used if you're trying to quickly work out if you have the right parameters for a particular dataset |


#### Using an interactive shell

Running `freyr` in an interactive shell is only recommended to check that parameters are valid or to run a test dataset, since your terminal will need to remain open and connected to BASC for the duration of the run, which could be many hours on real data.

To run the pipeline in an interactive shell, it is recommended to request a decent amount of time: at least 4-8 hours, perhaps even up to 24 hours, depending on the size of your dataset. Only a single CPU is needed for the `nextflow` 'head' job. A typical `sinteractive` command is below:

```
# this requests 1 CPU, 8 hours (480 minutes) of time, and 4 GB of memory
sinteractive -c 1 -t 480 --mem 4G
```

Once you have been allocated a computational node in the shell, run the commands as per the ['_Main pipeline command_'](/docs/insect_coi.md#main-pipeline-command) section above.

#### Using a SLURM script

SLURM scripts are the best way to run `freyr` for proper, large datasets, as your analysis will run on the server without needing interaction with your personal computer. BASC has a [good guide](http://users.basc.science.depi.vic.gov.au/jobs/slurm/sbatch/) to using `sbatch`, including a very useful [SLURM script generator](http://jsgen.basc.science.depi.vic.gov.au/). 

There is also a template SLURM script for running `freyr` on BASC available at [`./supplementary_scripts/basc_template.slurm`](/supplementary_scripts/basc/basc_template.slurm). 


### Optional: Running a test dataset

`freyr` comes packaged with a small (Illumina) test dataset that can be used to quickly check that the pipeline is working as expected and get you familiar with pipeline outputs. Typically, this should take only 5-10 minutes if the BASC job queue isn't full, and is most easily done in an interactive shell. To run `freyr` on this dataset, use the following commands:

```
## run these inside your cloned directory (ie. $working_dir)

# go to directory
cd $working_dir

# load Java
module load Java/17.0.6

# run nextflow on test data
NXF_VER=23.04.5 \
    nextflow run . \
    -profile basc_slurm,test

# once the dataset has run, clean up your analysis directory
rm -rf ./output/modules/* ./output/*.html ./work/*
```

### Monitoring a run of the pipeline

While `freyr` is running, something like the following will be displayed in the terminal (if in an interactive shell) or in the SLURM job output file (typically with a name like `slurm-*.out`):

```
Nextflow 24.04.4 is available - Please consider updating your version to it
N E X T F L O W  ~  version 23.04.5
Launching `./main.nf` [boring_hypatia] DSL2 - revision: 226d236840


:::::::::: :::::::::  :::::::::: :::   ::: :::::::::
:+:        :+:    :+: :+:        :+:   :+: :+:    :+:
+:+        +:+    +:+ +:+         +:+ +:+  +:+    +:+
:#::+::#   +#++:++#:  +#++:++#     +#++:   +#++:++#:
+#+        +#+    +#+ +#+           +#+    +#+    +#+
#+#        #+#    #+# #+#           #+#    #+#    #+#
###        ###    ### ##########    ###    ###    ###


~~~ freyr: Metabarcoding analysis for biosecurity and biosurveillance ~~~
 
Core Nextflow options
  runName        : boring_hypatia
  containerEngine: shifter
  launchDir      : /group/pathogens/IAWS/Personal/JackS/dev/freyr
  workDir        : /group/pathogens/IAWS/Personal/JackS/dev/freyr/work
  projectDir     : /group/pathogens/IAWS/Personal/JackS/dev/freyr
  userName       : js7t
  profile        : basc_slurm,test
  configFiles    : 

Main arguments
  samplesheet    : test_data/dual/samplesheet_read_dir.csv
  loci_params    : test_data/dual/loci_params.csv
  slurm_account  : pathogens

Max job request options
  max_cpus       : 1
  max_memory     : 2.GB
  max_time       : 10.m

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  slurm (1)
[21/225c1a] process > FREYR:PARSE_INPUTS (Whole dataset) [  0%] 0 of 1
[-        ] process > FREYR:SEQ_QC                       -
[-        ] process > FREYR:SPLIT_LOCI                   -
[-        ] process > FREYR:PRIMER_TRIM                  -
[-        ] process > FREYR:READ_FILTER                  -
[-        ] process > FREYR:FILTER_QUALPLOTS_PRE         -
[-        ] process > FREYR:FILTER_QUALPLOTS_POST        -
[-        ] process > FREYR:ERROR_MODEL_F                -
[-        ] process > FREYR:ERROR_MODEL_R                -
[-        ] process > FREYR:DENOISE1_F                   -
[-        ] process > FREYR:DENOISE1_R                   -
[-        ] process > FREYR:DADA_PRIORS_F                -
[-        ] process > FREYR:DADA_PRIORS_R                -
[-        ] process > FREYR:DENOISE2_F                   -
[-        ] process > FREYR:DENOISE2_R                   -
[-        ] process > FREYR:DADA_MERGEREADS              -
[-        ] process > FREYR:FILTER_SEQTAB                -
[-        ] process > FREYR:TAX_IDTAXA                   -
[-        ] process > FREYR:TAX_BLAST                    -
[-        ] process > FREYR:JOINT_TAX                    -
[-        ] process > FREYR:MERGE_TAX                    -
[-        ] process > FREYR:ASSIGNMENT_PLOT              -
[-        ] process > FREYR:TAX_SUMMARY                  -
[-        ] process > FREYR:TAX_SUMMARY_MERGE            -
[-        ] process > FREYR:PHYLOSEQ_UNFILTERED          -
[-        ] process > FREYR:PHYLOSEQ_FILTER              -
[-        ] process > FREYR:PHYLOSEQ_MERGE               -
[-        ] process > FREYR:READ_TRACKING                -

```

The bottom part of this output shows the status of each step of the pipeline (called a 'process' in `nextflow`) and can be used to judge how close to finishing the pipeline run is. 

If the pipeline completes without error, you should see something like this:

```
executor >  slurm (200)
[8d/6936dd] process > FREYR:PARSE_INPUTS (Whole dataset) [100%] 1 of 1 ✔
[43/879167] process > FREYR:SEQ_QC (K77JP)               [100%] 2 of 2 ✔
[69/d0540f] process > FREYR:SPLIT_LOCI (fwhF2-fwhR2nD... [100%] 16 of 16 ✔
[e1/505a7d] process > FREYR:PRIMER_TRIM (fwhF2-fwhR2n... [100%] 16 of 16 ✔
[f9/162c9b] process > FREYR:READ_FILTER (fwhF2-fwhR2n... [100%] 16 of 16 ✔
[5b/471d04] process > FREYR:FILTER_QUALPLOTS_PRE (K77... [100%] 16 of 16 ✔
[44/bb8661] process > FREYR:FILTER_QUALPLOTS_POST (K7... [100%] 16 of 16 ✔
[d8/bf502d] process > FREYR:ERROR_MODEL_F (EIF3LminiF... [100%] 4 of 4 ✔
[33/a48ec1] process > FREYR:ERROR_MODEL_R (fwhF2-fwhR... [100%] 4 of 4 ✔
[f2/8fe2d3] process > FREYR:DENOISE1_F (EIF3LminiF4-E... [100%] 16 of 16 ✔
[cb/720064] process > FREYR:DENOISE1_R (fwhF2-fwhR2nD... [100%] 16 of 16 ✔
[12/38a6a8] process > FREYR:DADA_PRIORS_F (EIF3LminiF... [100%] 4 of 4 ✔
[b1/bd1abc] process > FREYR:DADA_PRIORS_R (fwhF2-fwhR... [100%] 4 of 4 ✔
[24/39d9df] process > FREYR:DENOISE2_F (EIF3LminiF4-E... [100%] 16 of 16 ✔
[22/0d3023] process > FREYR:DENOISE2_R (fwhF2-fwhR2nD... [100%] 16 of 16 ✔
[0b/d2739f] process > FREYR:DADA_MERGEREADS (fwhF2-fw... [100%] 4 of 4 ✔
[60/6b60e8] process > FREYR:FILTER_SEQTAB (fwhF2-fwhR... [100%] 4 of 4 ✔
[8e/f6aca9] process > FREYR:TAX_IDTAXA (fwhF2-fwhR2nD... [100%] 4 of 4 ✔
[49/1db6d7] process > FREYR:TAX_BLAST (fwhF2-fwhR2nDa... [100%] 4 of 4 ✔
[0c/9280f5] process > FREYR:JOINT_TAX (fwhF2-fwhR2nDa... [100%] 4 of 4 ✔
[e0/99d8c1] process > FREYR:MERGE_TAX (EIF3LminiF4-EI... [100%] 2 of 2 ✔
[7c/830fbd] process > FREYR:ASSIGNMENT_PLOT (fwhF2-fw... [100%] 4 of 4 ✔
[7f/8d1b8b] process > FREYR:TAX_SUMMARY (fwhF2-fwhR2n... [100%] 4 of 4 ✔
[5d/5493c9] process > FREYR:TAX_SUMMARY_MERGE (Whole ... [100%] 1 of 1 ✔
[78/4f02d3] process > FREYR:PHYLOSEQ_UNFILTERED (EIF3... [100%] 2 of 2 ✔
[2b/be4016] process > FREYR:PHYLOSEQ_FILTER (EIF3Lmin... [100%] 2 of 2 ✔
[8b/a2b001] process > FREYR:PHYLOSEQ_MERGE (Whole dat... [100%] 1 of 1 ✔
[34/4c5996] process > FREYR:READ_TRACKING (Whole data... [100%] 1 of 1 ✔
[2024-08-08T14:57:07.284919555+10:00] >> Pipeline finished SUCCESSFULLY after 5m 31s
Completed at: 08-Aug-2024 14:57:08
Duration    : 5m 32s
CPU hours   : 0.6
Succeeded   : 200
```

At the end of the run, the pipeline will also produce two files in the `./output` directory: `report.html` and `timeline.html`. These detail the overall resource use of the pipeline and the time it took to run each process, respectively. If your run ends with an error, the `report.html` file will detail the error message at the top.

### Pipeline outputs

Currently, the outputs of the pipeline are found in the `./output/modules` directory, organised by process name. Most of these are not relevant for the typical user; useful final outputs are listed and explained below. 

#### Sequences, abundances, and taxonomic assignment

Inferred sequences (confusingly currently called both 'ASVs' and 'OTUs' in the pipeline) are filtered per locus/primer pair using the parameters specified in the `loci_params` file. As such, there are both filtered and unfiltered output files, which are found in the `./output/modules/phyloseq_unfiltered` and `./output/modules/phyloseq_filter` directories, respectively. Files in these directories are separated by locus, which is useful if your PCR primers target different genes. In contrast, the `./output/modules/phyloseq_merge` directory contains filtered and unfiltered results merged across primer pairs, which is useful if all your primer pairs target the same gene (eg. tagging technical or biological replicates). 

All three of these directories contain the following types of files (where `*` is a variable part of the file name):
- `asvs_*.fasta`: FastA file of OTU/ASV sequences
- `taxtab_*.csv`: taxonomic assignment of each OTU/ASV (name, no sequence) for every taxonomic rank
- `summary_*.csv`: abundance of OTUs/ASVs per sample (wide format), including sequence and taxonomic assignment
- `seqtab_*.csv`: simple table of OTU/ASV abundance per sample, not containing sequence or taxonomic assignment
- `samdf_*.csv`: samplesheets, split by primer pair, containing input metadata; samples removed during filtering will be missing from the `filtered` versions of this file
- `raw_*.csv`: large file of sample-level abundance per OTU/ASV (long format), including samplesheet metadata, loci parameters and taxonomic assignment (but not sequence)
- `ps_*.rds`: R data files in the [`phyloseq`](https://joey711.github.io/phyloseq/) format 

Detailed taxonomic assignment information for the unfiltered OTU/ASV set, merged across primer pairs, can also be found in `./output/modules/tax_summary_merge` in the `taxonomic_assignment_summary.csv` file. This includes the OTU/ASV name, sequence, primer pair, taxonomic assignment confidence per rank, and the identity of the top BLAST hit in the reference database. 

#### Visualisation and QC plots

The `phyloseq_unfiltered` directory contains **accumulation curve plots** (`accumulation_curve_*.pdf`) that show OTU/ASV accumulation curves per sample and flowcell. 

**Taxonomic assignment plots** (`*_taxonomic_assignment_summary.pdf`) can be found in `./output/modules/assignment_plot` and show the relationship between IDTAXA taxonomic assignment level and % identity to the top BLAST hit in the database, per flowcell and PCR primer combination. 

**OTU/ASV filtering plots** (`*_ASV_cleanup_summary.pdf`) can be found in `./output/modules/filter_seqtab` and show the number and abundance of sequences kept or filtered, within each flowcell and PCR primer combination.

A **read tracking plot** (`read_tracker.pdf`) can be found in `./output/modules/read_tracking`. This shows the total number of reads, per flowcell, that make it through each step of the pipeline, including the final taxonomic and abundance filtering of OTUs/ASVs. 

**Read quality plots** per sample and primer pair can be found in `./output/modules/filter_qualplots`, that visualise quality before (`*_pre_qualplots.pdf`) and after (`*_post_qualplots.pdf`) read filtering and truncation. These show the distribution of base quality per position (top) and the cumulative number of expected errors per read for various quantiles (bottom). A Nanopore-specific input QC report per sample from [`NanoPlot`](https://github.com/wdecoster/NanoPlot) (`*_NanoPlot-report.html`) can be found in `./output/modules/nanoplot`.

