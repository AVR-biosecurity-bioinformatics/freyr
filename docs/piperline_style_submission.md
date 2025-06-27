---
title: "Piperline style freyr submission"
output: html_document
date: "2025-01-22"
---

# Introduction

The below code is written for the Agriculture Victoria BASC computing cluster.

This workflow gives examples of how to run freyr through SLURM with a piperline style submission syntax


##Install nextflow in your home directory

freyr currently requires a specific version of nextflow that is not available as a BASC module. It also uses a third-party validation plugin that doesn't work with standalone (ie. module-based) distributions of nextflow. Luckily, it is very easy to install nextflow for a specific user on BASC:

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

NOTE: This only needs to be done once before any particular user uses nextflow for the first time -- you don't need to repeat this step for subsequent runs of this pipeline, or any other nextflow pipeline.


# Clone the freyr github repository

```
# Define the directory you will be running the analysis in
working_dir=/group/your.account/IAWS/Personal/Alexp/metabarcoding/marine_surveillance #CHANGE TO YOUR DIRECTORY

git clone https://github.com/AVR-biosecurity-bioinformatics/freyr $working_dir

```

## Updating freyr

To update to the latest version of freyr, run the below code in the terminal.

```
cd $working_dir
git pull

# if an error occurs, run:
git stash
git pull
```

# Demultiplex MiSeq run

For this workflow to run, we will need some sequencing runs to work with. If you are working with MiSeq data, it is recommended that the data is demultiplexed again using bcl2fastq, as the miseq does not put indexes in fasta headers by default which is required for the index swtiching calculation.

If you are not working with miseq data, the pipeline can still be run but the index switching calculations will be skipped

The below code is written for the Agriculture Victoria BASC computing cluster, and the locations will be different if you are using a different HPC cluster.

```
#load module
module load bcl2fastq2/2.20.0-GCC-13.3.0

#raise amount of available file handles
ulimit -n 4000

###Run1

#Set up input and outputs
inputdir=/group/sequencing/210219_M03633_0489_000000000-JDYG3 #CHANGE TO YOUR SEQ RUN

fcid=$(echo $inputdir | sed 's/^.*-//')
outputdir=$working_dir/data/$fcid
samplesheet=$inputdir/SampleSheet.csv

# convert samplesheet to unix format
dos2unix $samplesheet

#Demultiplex
bcl2fastq -p 12 --runfolder-dir $inputdir \
--output-dir $outputdir \
--sample-sheet $samplesheet \
--no-lane-splitting --barcode-mismatches 1

# Copy other necessary files and move fastqs
cd $outputdir
cp -r $inputdir/InterOp $outputdir
cp $inputdir/RunInfo.xml $outputdir
cp $inputdir/[Rr]unParameters.xml $outputdir
cp $samplesheet $outputdir
mv **/*.fastq.gz $outputdir

# Append fcid to start of sample names if missing
for i in *.fastq.gz; do
  if ! [[ $i == $fcid* ]]; then
  new=$(echo ${fcid} ${i}) #append together
  new=$(echo ${new// /_}) #remove any white space
  mv -v "$i" "$new"
  fi
done

```

## Copy reference databases

Reference databases for metaabrcoding are stored in /group/referencedata/mspd-db/metabarcoding/ on BASC

The latest terrestrial arthropod database is available here:
```
# Change 'folder-name' to the directory you are running the analysis in
cp /group/referencedata/mspd-db/metabarcoding/arthropod/terrestrial_arthropod_25_03_07/* $working_dir/reference/.

# Copy the PHMM model
cp /group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/folmer_fullength_model.rds $working_dir/reference/.
```

The older 'imappests' database is available here:
```
# Change 'folder-name' to the directory you are running the analysis in
cp /group/referencedata/mspd-db/metabarcoding/arthropod/imappests_coi_18_08_2020/* $working_dir/reference/.
```


# Submitting jobs


## Create primer parameters file

See [Usage](https://github.com/AVR-biosecurity-bioinformatics/freyr/blob/72-parameter-channels/docs/usage.md) for how to create a primer parameters file for your target loci


## Example 1: Single primer set

This is the 'classic' pipeRline style submission, where all fastq files, miseq sample sheets, and runparameters files are found automatically. 

This requires the data to be placed within the subfolders of the data directory. For example, if you are wishing to analyse two flow cells at once the data directory should look something like this:

    root/
      ├── data/
         ├── K739J/
         ├── JDYG3/
         
Paths to files can be relative to the submission location, or absolute paths     

**You will need to change your email address and account code below, see [BASC account code lists](https://users.basc.science.depi.vic.gov.au/jobs/slurm/slurm_accounts/)**

```
# Submit slurm job - CHANGE EMAIL TO YOUR OWN EMAIL
sbatch --mail-user=your.name@email.com --account=your.account supplementary_scripts/submit_slurm.sh \
--primer_params primer_params.csv \
--primers fwhF2-fwhR2n
```


# Example 2: Manually point to data folders

Alternatively, you can manually direct the submission script to the locations of oyur reads, sample sheets, and run parameters files. Multiple files or locations can be separated by a comma
This simple example 

if you wish to use multiple values (i.e. multiple reference databases, multiple primers per index) encapsulate them in a '' and separate them with a ;

Note, watch for the naming of the runParameters file, different versions of the miseq software output this file with a capitalized or lowercase R

**You will need to change your email address and account code below, see [BASC account code lists](https://users.basc.science.depi.vic.gov.au/jobs/slurm/slurm_accounts/)**

```
# Submit slurm job - CHANGE EMAIL TO YOUR OWN EMAIL
sbatch --mail-user=your.name@email.com --account=your.account supplementary_scripts/submit_slurm.sh \
--samplesheet 'data/K739J/SampleSheet_K739J.csv,data/JDYG3/SampleSheet_JDYG3.csv'  \
--primer_params primer_params.csv \
--read_dir 'data/K739J/,data/JDYG3/' \
--primers fwhF2-fwhR2n


# Note: to resume a previous run, add --nextflow-resume

```


**To resume a failed or paused run, add `--nextflow-resume` to the end of the command**


For more complex setups or for external users, please see [Usage](https://github.com/AVR-biosecurity-bioinformatics/freyr/blob/72-parameter-channels/docs/usage.md)

