#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads_paths, one or two read paths separated by ;
# $2 = sample_id
# $3 = seq_type (illumina, nanopore, pacbio)
# $4 = paired (true, false)

### define variables with better names

READS_PATHS=$1
SAMPLE_ID=$2
SEQ_TYPE=$3
PAIRED=$4


### process variables

if [ $SEQ_TYPE != "nanopore" ]; then
    echo "NANOPLOT module requires Nanopore reads -- check pipeline inputs and 'params.seq_type' parameter"
    exit 1
fi


# convert READS_PATHS into an array, conditional on presence of ; delimiter between paths
if [[ "$READS_PATHS" == *\;* ]]; then
    # IFS=';' read -r -a READ_ARRAY <<< "$READS_PATHS"
    # FWD_READS="${READ_ARRAY[0]}"
    # REV_READS="${READ_ARRAY[1]}"
    echo "Paired-end reads cannot be used with nanoplot; check 'params.paired' pipeline parameter"
    exit 1
else 
    SINGLE_READS=$READS_PATHS
fi

### run fastqc, with conditional code based on single vs paired

if [ $PAIRED == "true" ]; then 

    echo "Paired-end reads cannot be used with nanoplot; check 'params.paired' pipeline parameter"
    exit 1

    
elif [ $PAIRED == "false"  ] && [ $SEQ_TYPE == "nanopore" ]; then

    ### single-end reports
    NanoPlot \
        --fastq $SINGLE_READS \
        --tsv_stats \
        --prefix $SAMPLE_ID \
        -o .

else 
    echo "PAIRED variable must be true or false"
    exit 1
fi 

