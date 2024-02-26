#!/bin/bash
# Created by the DJPR job script generator for Slurm
# Mon Feb 12 2024 11:37:11 GMT+1100 (Australian Eastern Daylight Time)

# The name of the job:
#SBATCH --job-name="piperline.container_test"

# The project ID which this job should run under:
#SBATCH --account="pathogens"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The total amount of memory in megabytes in the job:
#SBATCH --mem=20480

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=jack.scanlan@agriculture.vic.gov.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi



# Run the job from this directory:
cd /group/pathogens/IAWS/Personal/JackS/piperline_tests

module purge

module load git/2.41.0-GCCcore-12.3.0-nodocs

git clone https://github.com/jackscanlan/piperline.git container_test && \
    cd container_test

module purge

cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/JDYG3 \
    ./data/

# replace old samplesheet with new samplesheet (just for this example)
rm ./data/JDYG3/SampleSheet.csv

cp /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/data/SampleSheet_JDYG3.csv \
    ./data/JDYG3/SampleSheet.csv

cp -r /group/pathogens/IAWS/Projects/Metabarcoding/dros_surveillance/reference \
    .

module load shifter

# shifterimg pull docker:jackscanlan/piperline:0.0.2

shifter --image jackscanlan/piperline:0.0.2 \
    Rscript ./jack_notes/basc_shifter.R

exit