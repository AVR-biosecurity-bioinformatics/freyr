#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads_paths, one or two read paths separated by ;
# $2 = seq_type (illumina, nanopore, pacbio)
# $3 = paired (true, false)

### define variables with better names

READS_PATHS=$1
SEQ_TYPE=$2
PAIRED=$3


### process variables

# convert READS_PATHS into an array, conditional on presence of ; delimiter between paths
if [[ "$READS_PATHS" == *\;* ]]; then
    IFS=';' read -r -a READ_ARRAY <<< "$READS_PATHS"
    FWD_READS="${READ_ARRAY[0]}"
    REV_READS="${READ_ARRAY[1]}"
else 
    SINGLE_READS=$READS_PATHS
fi

### run fastqc, with conditional code based on single vs paired

if [ $PAIRED == "true" ]; then 

    ### paired-end reports
    fastqc \
    -t 1 \
    --outdir . \
    $FWD_READS \
    $REV_READS

    
elif [ $PAIRED == "false" ]; then

    ### single-end reports
        fastqc \
    -t 1 \
    --outdir . \
    $SINGLE_READS

else 
    echo "PAIRED variable must be true or false"
    exit 1
fi 

