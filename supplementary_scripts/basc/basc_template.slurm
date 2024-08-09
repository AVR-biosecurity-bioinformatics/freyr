#!/bin/bash

# The name of the job:
#SBATCH --job-name="freyr"
# you can change this to whatever you like

# The name of the output file (name of the job name plus the job ID):
#SBATCH --output="%x_%j.out"

# The project ID which this job should run under:
#SBATCH --account="pathogens"
# change to a different account if you're not part of the pathogens group

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The total amount of memory for the job:
#SBATCH --mem=4G

# Send yourself an email when the job aborts, begins and ends:
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=<your-email-address>@agriculture.vic.gov.au
# change to your internal email address

# The maximum running time of the job in days-hours:mins:sec (this runs for 12 hours)
#SBATCH --time=0-12:0:00
# change to longer if you have a very large dataset (ie. >2-3 flowcells, very deep NovaSeq sequencing)

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

#### running freyr pipeline

## NOTE: run a new, modified version of this SLURM script after you have cloned the freyr GitHub repository into a new working directory,
## and created and uploaded your samplesheet and loci_params .csv files to the working directory.
## For a guide to these steps, see here: https://github.com/AVR-biosecurity-bioinformatics/freyr/blob/main/docs/insect_coi.md

# define your working (analysis) directory; recommended to be an absolute path (ie. starts with '/group/')
WORKING_DIR=/path/to/your/working/directory # make sure you change this

# define your samplesheet and loci_params file paths (can be absolute, or relative to $WORKING_DIR)
SAMPLESHEET=/path/to/samplesheet.csv
LOCI_PARAMS=/path/to/loci_params.csv

# change to working directory
cd $WORKING_DIR

# remove old output and work files if you have done previous runs 
rm -rf ./output/modules/* ./output/*.html ./work/*

# load Java module
module load Java/17.0.6

# run pipeline with nextflow
NXF_VER=23.04.5 \
    nextflow run . \
    --samplesheet $SAMPLESHEET \
    --loci_params $LOCI_PARAMS \
    -profile basc_slurm
# add in any additional pipeline parameters (eg. --miseq_internal) after the last line above
# don't forget to escape new lines with '\' like above

# exit script
exit


#### submitting this script as a job

## save this script in $WORKING_DIR/supplementary_scripts with a new filename, then run the below command (uncommented with new filename) while in $WORKING_DIR:
# sbatch $WORKING_DIR/supplementary_scripts/<this-file>.slurm