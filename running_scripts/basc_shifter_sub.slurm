#!/bin/bash

# The name of the job:
#SBATCH --job-name="piperline.container_test"

# The name of the output file (name of the job name plus the job ID):
#SBATCH --output="%x_%j.out"

# The project ID which this job should run under:
#SBATCH --account="pathogens"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The total amount of memory in megabytes in the job:
#SBATCH --mem=20480

# Send yourself an email when the job aborts, begins and ends:
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=jack.scanlan@agriculture.vic.gov.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

#### Location of run directory ####
RUN_LOC=/group/pathogens/IAWS/Personal/JackS/piperline_tests

#### Name of run directory; to be made ####
RUN_DIR=subsample_test

# move to the specified directory after checking that directory exists and is not a file
if [ -d "${RUN_LOC:+$RUN_LOC/}" ]; then
    cd $RUN_LOC
else
    "Run directory is not a valid directory"
fi 

# clone GitHub repo for the pipeline and move into it
module purge
module load git/2.41.0-GCCcore-12.3.0-nodocs
git clone https://github.com/jackscanlan/piperline.git $RUN_DIR && \
    cd $RUN_DIR
module purge

# copy subsample data into running directory
cp -r /group/pathogens/IAWS/Personal/JackS/piperline_tests/test_data/subsamples \
    ./data/

# # replace old samplesheet with new samplesheet
# rm ./data/JDYG3/SampleSheet.csv
# cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv \
#     ./data/JDYG3/SampleSheet.csv

# copy across existing reference database
cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference \
    .

# run pipeline using Docker container
module load shifter
# shifterimg pull docker:jackscanlan/piperline:0.0.2 ## uncomment if needing to pull image
shifter --image jackscanlan/piperline:0.0.2
Rscript ./running_scripts/basc_shifter.R

# exit script
exit