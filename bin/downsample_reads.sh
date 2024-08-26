#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads_paths, one or two read paths separated by ;
# $2 = seq_type (illumina, nanopore, pacbio)
# $3 = paired (true, false)
# $4 = params.downsample_reads (integer)

### define variables with better names

READS_PATHS=$1
SEQ_TYPE=$2
PAIRED=$3
DOWNSAMPLE_READS=$4


### process variables

SEED=1

if [ $DOWNSAMPLE_READS == "null" ]; then
    echo "DOWNSAMPLE_READS module requires 'params.downsample_reads' to be an integer >=1 -- check pipeline inputs and 'params.downsample_reads' parameter"
    exit 1
fi


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
    
    ### paired-end downsampling
    ## forward reads
    seqtk sample \
        -s $SEED \
        $FWD_READS \
        $DOWNSAMPLE_READS \
        > "down_${FWD_READS}"

    ## reverse reads
    seqtk sample \
        -s $SEED \
        $REV_READS \
        $DOWNSAMPLE_READS \
        > "down_${REV_READS}"

    
elif [ $PAIRED == "false"  ] && [ $SEQ_TYPE == "nanopore" ]; then

    ### single-end downsampling
    seqtk sample \
        -s $SEED \
        $SINGLE_READS \
        $DOWNSAMPLE_READS \
        > "down_${SINGLE_READS}"

else 
    echo "PAIRED variable must be true or false"
    exit 1
fi 

