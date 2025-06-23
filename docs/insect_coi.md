# Insect COI workflow for AgVic researchers

#### by Jack Scanlan, updated 2025-06**

This tutorial will help you run `freyr` on the Agriculture Victoria BASC HPC system, for experiments utilising the 'freshwater' COI primers `fwhF2-fwhR2n` on mixed arthropod samples (focusing on insects).

## Prior to analysis

The following are some assumptions that you need to make sure are met before you start this analysis:

1. Samples have been sequenced with Illumina paired-end sequencing (eg. MiSeq or NovaSeq)
2. You have used the `fwhF2-fwhR2n` primers to amplify the COI gene, with Illumina adapters added to the end of each amplicon, ie. amplicons have **not** been fragmented before adapters were added
3. Your data is available, per sample, as gzipped FastQ files (with `.fastq.gz` or `.fq.gz` extensions), with forward and reverse reads in separate files, ie. **not** interleaved

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

| Field | Meaning | Requirements | Change to: |
| --- | --- | --- | --- |
| `sample` | Unique sample ID | Must be unique within the samplesheet. If `read_dir` is used over `fwd`/`rev`/`single`, it must be present at the start of the name of each read file. | Depends on data |
| `read_group` | Sequencing read group | If this is unknown, do not include this field and it will be automatically determined from the read headers of each sample. | Depends on data  |
| `primers` | Name of the PCR primer pair(s) | Can specify multiple by separating with a semi-colon. Must be a subset of the values in the `--primer_params` file. Each sample can have a different value in this field if appropriate. | `fwhF2-fwhR2n` |
| `read_dir` | Directory containing sequencing reads | If this is specified, the pipeline searches this directory (including subdirectories) for read files starting with `sample`. Cannot be used in conjunction with `fwd` and `rev`; delete field if not using. | Up to you! |
| `fwd` | Exact path of the forward read file for this sample | Cannot be used in conjunction with `read_dir`; delete field if not using. | Up to you!  |
| `rev` | Exact path of the reverse read file for this sample | Cannot be used in conjunction with `read_dir`; delete field if not using. | Up to you!  |

The easiest way to specify the reads for each sample is to copy your read files into the `./data` directory and then put `./data` in the `read_dir` field for every row in the samplesheet. \

> IMPORTANT NOTE: You must only use one of either `read_dir` OR `fwd` + `rev` in the samplesheet. Delete the fields/columns you're not using.

> NOTE: Paths for the `read_dir` and `fwd`/`rev` fields can be absolute or relative.

Once you have made your samplesheet `.csv` (it can have any name), upload it to the `inputs` directory.

### Make your primer parameters file

The primer parameters (or `primer_params`) file is a `.csv` that conveys information to the pipeline about how to process amplicons derived from different PCR primers. One row should be used per primer pair. A 'default' file can be found at `./inputs/primer_params_default.csv` that contains defaults for 'typical' COI runs, although you should always make sure you don't need to change any of these.

Every field in the `primer_params` file is currently required and must be specified, so don't remove any of the columns.

For the purposes of this run, we only need to change a few fields in `primer_params_default.csv`:

| field | meaning | change to: |
| --- | --- | --- |
| `read_trunc_length` | Number of bp to truncate each read to | If your read length is 250 or 251 bp (ie. typical MiSeq run), keep as `150`; if read lengths are higher or lower, you might need to adjust accordingly so you both remove poor sequence quality at the end of reads **and** keep overlap between paired reads |
| `phmm` | Path to PHMM model of COI | `/group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/folmer_fullength_model.rds` |
| `idtaxa_db` | Path to trained IDTAXA database | `/group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/idtaxa_bftrimmed.rds` |
| `ref_fasta` | Path to reference sequence database  | `/group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/insecta_hierarchial_bftrimmed.fa.gz` |

Once you have made your primer parameters `.csv` (it can have any name), upload it to the `inputs` directory.

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
NXF_VER=23.05.0-edge \
    nextflow run . \
    --samplesheet ./inputs/[your_samplesheet].csv \
    --primer_params ./inputs/[your_primer_params].csv \
    -profile basc_slurm

```

The above code is doing the following:

- `NXF_VER=23.05.0-edge` tells `nextflow` to run as version 23.05.0-edge, which is currently required for the pipeline to work correctly on BASC
- `nextflow run .` tells `nextflow` to run the pipeline using files in the current directory
- `--samplesheet` and `--primer_params` specify where your samplesheet and primer parameters files are; these can be relative or absolute paths, here we use relative paths
- `-profile basc_slurm` tells `nextflow` to run using settings that make sense for BASC; in particular, it will spawn a new job for every process, making its computational execution very efficient, and will use `shifter` to run Docker containers for all the software tools

> NOTE: Confusingly, pipeline parameters like `--samplesheet` need to be typed on the command line using a double hyphen, while `nextflow` parameters like `-profile` need to be typed with a single hyphen. When in doubt, copy-paste from this page.

#### Advanced pipeline parameters

`freyr` has optional, advanced parameters that can be specified on the command line after `-profile`. Below are some relevant ones, but you typically won't need to alter these.

| parameter | meaning | requirements |
| --- | --- | --- |
| `--miseq_internal` | If your sequencing data was generated by a MiSeq machine and re-demultiplexed to put indices in headers, use this option to estimate index switching and some other QC stats | Can be `false` (default) or `true`; ignore if you're uncertain about this |
| `--primer_error_rate` | Error rate allowed when detecting primer sequences | Sets the value for `cutadapt -e`; see the [`cutadapt` reference manual](https://cutadapt.readthedocs.io/en/v4.7/guide.html#error-tolerance) for more info |
| `--primer_n_trim` | Recognise N bases in reads when detecting primers | Can be `false` (default) or `true`; only use this if your sequencing data quality is very low |
| `--subsample` | Reduce the number of input samples per `primers` x `read_group` combination to this number | Should be a number >=1; this is mainly used for development and debugging and should only be used if you're trying to quickly work out if you have the right parameters for a particular dataset |

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

`freyr` comes packaged with a small test dataset that can be used to quickly check that the pipeline is working as expected and get you familiar with pipeline outputs. Typically, this should take only 5-10 minutes if the BASC job queue isn't full, and is most easily done in an interactive shell. To run `freyr` on this dataset, use the following commands:

```
## run these inside your cloned directory (ie. $working_dir)

# go to directory
cd $working_dir

# load Java
module load Java/17.0.6

# run nextflow on test data
NXF_VER=23.05.0-edge \
    nextflow run . \
    -profile basc_slurm,test

# once the dataset has run, clean up your analysis directory
rm -rf ./output/modules/* ./output/run_info/* ./output/results/* ./work/*
```

### Monitoring a run of the pipeline

While `freyr` is running, something like the following will be displayed in the terminal (if in an interactive shell) or in the SLURM job output file (typically with a name like `slurm-*.out`):

```
Nextflow 25.05.0-edge is available - Please consider updating your version to it
N E X T F L O W  ~  version 23.05.0-edge
Launching `./main.nf` [goofy_rubens] DSL2 - revision: 690caecd45


:::::::::: :::::::::  :::::::::: :::   ::: :::::::::
:+:        :+:    :+: :+:        :+:   :+: :+:    :+:
+:+        +:+    +:+ +:+         +:+ +:+  +:+    +:+
:#::+::#   +#++:++#:  +#++:++#     +#++:   +#++:++#:
+#+        +#+    +#+ +#+           +#+    +#+    +#+
#+#        #+#    #+# #+#           #+#    #+#    #+#
###        ###    ### ##########    ###    ###    ###


~~~ freyr: Metabarcoding analysis for biosecurity and biosurveillance ~~~
 
Core Nextflow options
  runName        : goofy_rubens
  containerEngine: shifter
  launchDir      : /group/pathogens/IAWS/Personal/JackS/dev/freyr
  workDir        : /group/pathogens/IAWS/Personal/JackS/dev/freyr/work
  projectDir     : /group/pathogens/IAWS/Personal/JackS/dev/freyr
  userName       : js7t
  profile        : basc_slurm,test,debug
  configFiles    : 

Main arguments
  samplesheet    : test_data/dual/samplesheet_read_dir.csv
  primer_params  : test_data/dual/primer_params.csv
  miseq_internal : true
  miseq_dir      : test_data/dual
  slurm_account  : ngdsi

Debugging options
  rdata          : true

Max job request options
  max_cpus       : 1
  max_memory     : 2.GB
  max_time       : 10.m

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
executor >  slurm (1)
[6d/1a415a] process > FREYR:PARSE_INPUTS (Whole dataset)         [  0%] 0 of 1
[-        ] process > FREYR:PROCESS_READS:MISEQ_QC               -
[-        ] process > FREYR:PROCESS_READS:FASTQC                 -
[-        ] process > FREYR:PROCESS_READS:SPLIT_LOCI             -
[-        ] process > FREYR:PROCESS_READS:PRIMER_TRIM            -
[-        ] process > FREYR:PROCESS_READS:READ_FILTER            -
[-        ] process > FREYR:PROCESS_READS:FILTER_QUALPLOTS_PRE   -
[-        ] process > FREYR:PROCESS_READS:FILTER_QUALPLOTS_POST  -
[-        ] process > FREYR:DADA2:ERROR_MODEL_F                  -
[-        ] process > FREYR:DADA2:ERROR_MODEL_R                  -
[-        ] process > FREYR:DADA2:DENOISE1_F                     -
[-        ] process > FREYR:DADA2:DENOISE1_R                     -
[-        ] process > FREYR:DADA2:PRIORS_F                       -
[-        ] process > FREYR:DADA2:DENOISE2_F                     -
[-        ] process > FREYR:DADA2:PRIORS_R                       -
[-        ] process > FREYR:DADA2:DENOISE2_R                     -
[-        ] process > FREYR:DADA2:MAKE_SEQTAB_PAIRED             -
[-        ] process > FREYR:FILTERING:FILTER_CHIMERA             -
[-        ] process > FREYR:FILTERING:FILTER_LENGTH              -
[-        ] process > FREYR:FILTERING:FILTER_PHMM                -
[-        ] process > FREYR:FILTERING:FILTER_FRAME               -
[-        ] process > FREYR:FILTERING:MERGE_FILTERS              -
[-        ] process > FREYR:TAXONOMY:TAX_IDTAXA                  -
[-        ] process > FREYR:TAXONOMY:TAX_BLAST                   -
[-        ] process > FREYR:TAXONOMY:JOINT_TAX                   -
[-        ] process > FREYR:TAXONOMY:MERGE_TAX                   -
[-        ] process > FREYR:TAXONOMY:ASSIGNMENT_PLOT             -
[-        ] process > FREYR:TAXONOMY:TAX_SUMMARY                 -
[-        ] process > FREYR:TAXONOMY:TAX_SUMMARY_MERGE           -
[-        ] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_UNFILTERED -
[-        ] process > FREYR:RESULT_SUMMARIES:ACCUMULATION_CURVE  -
[-        ] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_FILTER     -
[-        ] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_MERGE      -
[-        ] process > FREYR:RESULT_SUMMARIES:READ_TRACKING       -

```

The bottom part of this output shows the status of each step of the pipeline (called a 'process' in `nextflow`) and can be used to judge how close to finishing the pipeline run is. 

If the pipeline completes without error, you should see something like this:

```
executor >  slurm (216)
[6d/1a415a] process > FREYR:PARSE_INPUTS (Whole dataset)                                              [100%] 1 of 1 ✔
[4b/1fa863] process > FREYR:PROCESS_READS:MISEQ_QC (000000000-K739J__1)                               [100%] 2 of 2 ✔
[0b/0a876c] process > FREYR:PROCESS_READS:FASTQC (K77JP_Mock2)                                        [100%] 8 of 8 ✔
[bd/28ead2] process > FREYR:PROCESS_READS:SPLIT_LOCI (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)            [100%] 16 of 16 ✔
[59/3b4973] process > FREYR:PROCESS_READS:PRIMER_TRIM (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)           [100%] 16 of 16 ✔
[7b/d61555] process > FREYR:PROCESS_READS:READ_FILTER (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)           [100%] 16 of 16 ✔
[2b/677c74] process > FREYR:PROCESS_READS:FILTER_QUALPLOTS_PRE (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)  [100%] 16 of 16 ✔
[7c/822e26] process > FREYR:PROCESS_READS:FILTER_QUALPLOTS_POST (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4) [100%] 16 of 16 ✔
[3c/b4ad4c] process > FREYR:DADA2:ERROR_MODEL_F (fwhF2-fwhR2nDac; 000000000-K77JP__1)                 [100%] 4 of 4 ✔
[88/a28e70] process > FREYR:DADA2:ERROR_MODEL_R (EIF3LminiF4-EIF3lminiR4; 000000000-K77JP__1)         [100%] 4 of 4 ✔
[9c/6b9305] process > FREYR:DADA2:DENOISE1_F (K77JP_Trap1_fwhF2-fwhR2nDac)                            [100%] 16 of 16 ✔
[27/1fa697] process > FREYR:DADA2:DENOISE1_R (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)                    [100%] 16 of 16 ✔
[ec/3689b3] process > FREYR:DADA2:PRIORS_F (fwhF2-fwhR2nDac; 000000000-K77JP__1)                      [100%] 4 of 4 ✔
[02/0a6407] process > FREYR:DADA2:DENOISE2_F (K77JP_Trap1_fwhF2-fwhR2nDac)                            [100%] 16 of 16 ✔
[f5/26ea6c] process > FREYR:DADA2:PRIORS_R (EIF3LminiF4-EIF3lminiR4; 000000000-K77JP__1)              [100%] 4 of 4 ✔
[78/6911ce] process > FREYR:DADA2:DENOISE2_R (K77JP_Trap2_EIF3LminiF4-EIF3lminiR4)                    [100%] 16 of 16 ✔
[18/3b4307] process > FREYR:DADA2:MAKE_SEQTAB_PAIRED (EIF3LminiF4-EIF3lminiR4; 000000000-K739J__1)    [100%] 4 of 4 ✔
[73/7f7edf] process > FREYR:FILTERING:FILTER_CHIMERA (EIF3LminiF4-EIF3lminiR4)                        [100%] 2 of 2 ✔
[db/3f8ade] process > FREYR:FILTERING:FILTER_LENGTH (EIF3LminiF4-EIF3lminiR4)                         [100%] 2 of 2 ✔
[79/4c188c] process > FREYR:FILTERING:FILTER_PHMM (EIF3LminiF4-EIF3lminiR4)                           [100%] 2 of 2 ✔
[df/7721df] process > FREYR:FILTERING:FILTER_FRAME (EIF3LminiF4-EIF3lminiR4)                          [100%] 2 of 2 ✔
[e4/7d4e66] process > FREYR:FILTERING:MERGE_FILTERS (EIF3LminiF4-EIF3lminiR4)                         [100%] 2 of 2 ✔
[8d/c2adc3] process > FREYR:TAXONOMY:TAX_IDTAXA (fwhF2-fwhR2nDac; 000000000-K77JP__1)                 [100%] 4 of 4 ✔
[0c/6ce163] process > FREYR:TAXONOMY:TAX_BLAST (fwhF2-fwhR2nDac; 000000000-K739J__1)                  [100%] 4 of 4 ✔
[d0/fae64a] process > FREYR:TAXONOMY:JOINT_TAX (EIF3LminiF4-EIF3lminiR4; 000000000-K739J__1)          [100%] 4 of 4 ✔
[59/710335] process > FREYR:TAXONOMY:MERGE_TAX (fwhF2-fwhR2nDac)                                      [100%] 2 of 2 ✔
[c8/495a60] process > FREYR:TAXONOMY:ASSIGNMENT_PLOT (EIF3LminiF4-EIF3lminiR4; 000000000-K739J__1)    [100%] 4 of 4 ✔
[8b/1b03e4] process > FREYR:TAXONOMY:TAX_SUMMARY (EIF3LminiF4-EIF3lminiR4; 000000000-K739J__1)        [100%] 4 of 4 ✔
[f7/2435f6] process > FREYR:TAXONOMY:TAX_SUMMARY_MERGE (Whole dataset)                                [100%] 1 of 1 ✔
[c4/db5c8e] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_UNFILTERED (fwhF2-fwhR2nDac)                    [100%] 2 of 2 ✔
[c6/4bb480] process > FREYR:RESULT_SUMMARIES:ACCUMULATION_CURVE (EIF3LminiF4-EIF3lminiR4)             [100%] 2 of 2 ✔
[a1/16f2fd] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_FILTER (fwhF2-fwhR2nDac)                        [100%] 2 of 2 ✔
[83/f2051f] process > FREYR:RESULT_SUMMARIES:PHYLOSEQ_MERGE (Whole dataset)                           [100%] 1 of 1 ✔
[0e/9e3059] process > FREYR:RESULT_SUMMARIES:READ_TRACKING (Whole dataset)                            [100%] 1 of 1 ✔
[2025-06-23T11:03:04.820843603+10:00] >> Pipeline finished SUCCESSFULLY after 20m 9s
Completed at: 23-Jun-2025 11:03:06
Duration    : 20m 10s
CPU hours   : 0.6
Succeeded   : 216
```

At the end of the run, the pipeline will also produce two files in the `./output` directory: `report.html` and `timeline.html`. These detail the overall resource use of the pipeline and the time it took to run each process, respectively. If your run ends with an error, the `report.html` file will detail the error message at the top.

### Pipeline outputs

Currently, the outputs of the pipeline are found in the `./output/modules` directory, organised by process name. Most of these are not relevant for the typical user; useful final outputs are listed and explained below.

#### Sequences, abundances, and taxonomic assignment

Inferred sequences (ASVs) are filtered per locus/primer pair using the parameters specified in the `primer_params` file. As such, there are both filtered and unfiltered output files, which are found in the `./output/modules/phyloseq_unfiltered` and `./output/modules/phyloseq_filter` directories, respectively. Files in these directories are separated by locus, which is useful if your PCR primers target different genes. In contrast, the `./output/modules/phyloseq_merge` directory contains filtered and unfiltered results merged across primer pairs, which is useful if all your primer pairs target the same gene (eg. tagging technical or biological replicates).

All three of these directories contain the following types of files (where `*` is a variable part of the file name):

- `asvs_*.fasta`: FASTA file of ASV sequences
- `taxtab_*.csv`: taxonomic assignment of each ASV (name, no sequence) for every taxonomic rank
- `summary_*.csv`: abundance of ASVs per sample (wide format), including sequence and taxonomic assignment
- `seqtab_*.csv`: simple table of ASV abundance per sample, not containing sequence or taxonomic assignment
- `samdf_*.csv`: samplesheets, split by primer pair, containing input metadata; samples removed during filtering will be missing from the `filtered` versions of this file
- `raw_*.csv`: large file of sample-level abundance per ASV (long format), including samplesheet metadata, primer parameters, taxonomic assignment and sequence
- `ps_*.rds`: R data files in the [`phyloseq`](https://joey711.github.io/phyloseq/) format

Detailed taxonomic assignment information for the unfiltered ASV set, merged across primer pairs, can also be found in `./output/modules/tax_summary_merge` in the `taxonomic_assignment_summary.csv` file. This includes the ASV name, sequence, primer pair, taxonomic assignment confidence per rank, and the identity of the top BLAST hit in the reference database.

#### Visualisation and QC plots

The `accumulation_curve` directory contains **accumulation curve plots** (`accumulation_curve_*.pdf`) that show ASV accumulation curves per sample for a particular primer pair. 

**Taxonomic assignment plots** (`*_taxonomic_assignment_summary.pdf`) can be found in `./output/modules/assignment_plot` and show the relationship between IDTAXA taxonomic assignment level and % identity to the top BLAST hit in the database, per flowcell and PCR primer combination.

**ASV filtering plots** (`*_ASV_cleanup_summary.pdf`) can be found in `./output/modules/filter_seqtab` and show the number and abundance of sequences kept or filtered, within each flowcell and PCR primer combination.

A **read tracking plot** (`read_tracker.pdf`) can be found in `./output/modules/read_tracking`. This shows the total number of reads, per flowcell, that make it through each step of the pipeline, including the final taxonomic and abundance filtering of ASVs.

**Read quality plots** per sample and primer pair can be found in `./output/modules/filter_qualplots`, that visualise quality before (`*_pre_qualplots.pdf`) and after (`*_post_qualplots.pdf`) read filtering and truncation. These show the distribution of base quality per position (top) and the cumulative number of expected errors per read for various quantiles (bottom).

If the pipeline parameter `--miseq_internal` was set to `true` and the data directory had the correct MiSeq-output structure, there will be **additional QC plots** in `./output/modules/seq_qc` per flowcell. `*_flowcell_qc.pdf` has stats about the quality of the MiSeq runs, while `*_index_switching.pdf` shows a heat map of the calculated index switching rate per index combination.