#!/bin/bash

# The name of the job:
#SBATCH --job-name="bs.full_teph"

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
cd /group/pathogens/IAWS/Personal/JackS/nextflow_tests

# clone github repo 
module purge
module load git/2.41.0-GCCcore-12.3.0-nodocs
git clone https://github.com/jackscanlan/piperline.git full_teph && \
    cd full_teph
module purge

# copy data
cp -r /group/pathogens/IAWS/Personal/JackS/nextflow_tests/piperline_nextflow/test_data/full_teph/* \
    ./data/

# copy reference dbs
cp -r /group/pathogens/IAWS/Personal/JackS/nextflow_tests/piperline_nextflow/reference \
    .

# run code in container
module load shifter
# shifterimg pull docker:jackscanlan/piperline:0.0.2
shifter --image jackscanlan/piperline:0.0.2 \
    Rscript ./jack_notes/bs_full_teph.R

exit