#!/bin/bash
# Created by the DJPR job script generator for Slurm
# Mon Feb 12 2024 11:37:11 GMT+1100 (Australian Eastern Daylight Time)

# The name of the job:
#SBATCH --job-name="piperline.dros_test"

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
#SBATCH --time=2-0:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from this directory:
cd /group/pathogens/IAWS/Personal/JackS/piperline_tests/dros_test


# This is required to ensure that modules loaded in your shell will not affect the job
# All modules required for the job need to be loaded by this script 
module purge

####

# load modules
module load R/4.2.0-foss-2021b
module load pkgconfig/1.5.1-GCCcore-9.3.0-Python-3.8.2
module load GDAL/3.3.0-foss-2021a
module load BLAST+/2.11.0-gompi-2020a
module load Pandoc/2.5
module load ZeroMQ/4.3.2-GCCcore-9.3.0

# run script that runs pipeline
Rscript ./basc_run.R

