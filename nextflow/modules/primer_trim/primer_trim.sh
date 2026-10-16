#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads_paths, one or two read paths separated by ;
# $2 = process_params.for_primer_seq, aka. fwd primer sequence
# $3 = process_params.rev_primer_seq, aka. rev primer sequence
# $4 = primers, aka. name of PCR primer pair
# $5 = process_params.locus, aka. name of target locus
# $6 = sample_primers, aka. sample ID
# $7 = read_group
# $8 = process_params.primer_n_trim (true or false)
# $9 = process_params.primer_error_rate
# $10 = process_params.seq_type (illumina, nanopore, pacbio)
# $11 = process_params.paired (true, false)

### define variables with better names

READS_PATHS=$1
FOR_PRIMER_SEQ=$2
REV_PRIMER_SEQ=$3
PRIMERS=$4
LOCUS=$5
SAMPLE_PRIMERS=$6
READ_GROUP=$7
PRIMER_N_TRIM=$8
PRIMER_ERROR_RATE=$9
SEQ_TYPE=${10}
PAIRED=${11}


### process variables

# convert READS_PATHS into an array, conditional on presence of ; delimiter between paths
if [[ "$READS_PATHS" == *\;* ]]; then
    IFS=';' read -r -a READ_ARRAY <<< "$READS_PATHS"
    FWD_READS="${READ_ARRAY[0]}"
    REV_READS="${READ_ARRAY[1]}"
else 
    SINGLE_READS=$READS_PATHS
fi

# change "I" in primer seq to "N", if present
FWD_PRIMER=${2/I/N}
REV_PRIMER=${3/I/N}

# reverse complement the primer sequences
FWD_PRIMER_RC=$(echo ${FWD_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)
REV_PRIMER_RC=$(echo ${REV_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)

### optional flags

# parse 'params.primer_n_trim' into cutadapt flag
if [ $PRIMER_N_TRIM == "true" ]; then
    PRIMER_N_TRIM_ARG="--match-read-wildcards"
else 
    PRIMER_N_TRIM_ARG=""
fi

# combine optional flags
OPTIONAL_ARGS=$( echo $PRIMER_N_TRIM_ARG )

### trim primers, with conditional code based on single vs paired

if [ $PAIRED == "true" ]; then 

    ### paired-end trimming

    ## run cutadapt in "linked adapter" mode
    # see here: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads
    cutadapt \
        -a ^${FWD_PRIMER}...${REV_PRIMER_RC} \
        -A ^${REV_PRIMER}...${FWD_PRIMER_RC} \
        --discard-untrimmed \
        --rename="{header}" \
        --report=minimal \
        --minimum-length 20 \
        --times 2 \
        --revcomp \
        -e $PRIMER_ERROR_RATE \
        $OPTIONAL_ARGS \
        -o ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R1.fastq.gz \
        -p ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R2.fastq.gz \
        ${FWD_READS} \
        ${REV_READS}

    ## count reads in output files
    # forward reads
    if [ -f ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R1.fastq.gz ]; then
        R1_OUT_LINES=$(zcat ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R1.fastq.gz | wc -l)
        R1_OUT=$(( $R1_OUT_LINES / 4 ))
    else 
        R1_OUT=0
    fi

    # reverse reads
    if [ -f ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R2.fastq.gz ]; then
        R2_OUT_LINES=$(zcat ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R2.fastq.gz | wc -l)
        R2_OUT=$(( $R2_OUT_LINES / 4 ))
    else 
        R2_OUT=0
    fi

    # save as .csv 
    echo "primer_trim,$SAMPLE_PRIMERS,$READ_GROUP,$PRIMERS,$R1_OUT,$R2_OUT" > primer_trim_${SAMPLE_PRIMERS}_${PRIMERS}_readsout.csv # columns: stage; sample_primers; read_group; primers; fwd_out; fwd_out

elif [ $PAIRED == "false" ]; then 
    
    ### single-end trimming

    ## run cutadapt in "linked adapter" mode
    # see here: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads
    cutadapt \
        -g X${FWD_PRIMER}...${REV_PRIMER_RC}X \
        --discard-untrimmed \
        --rename="{header}" \
        --report=minimal \
        --minimum-length 20 \
        --times 2 \
        --revcomp \
        -e $PRIMER_ERROR_RATE \
        $OPTIONAL_ARGS \
        -o ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R0.fastq.gz \
        ${SINGLE_READS} 
    
    ## count reads in output files
    # single reads
    if [ -f ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R0.fastq.gz ]; then
        R0_OUT_LINES=$(zcat ${SAMPLE_PRIMERS}_${LOCUS}_${PRIMERS}_trim_R0.fastq.gz | wc -l)
        R0_OUT=$(( $R0_OUT_LINES / 4 ))
    else 
        R0_OUT=0
    fi

    # save as .csv 
    echo "primer_trim,$SAMPLE_PRIMERS,$READ_GROUP,$PRIMERS,$R0_OUT,$R0_OUT" > primer_trim_${SAMPLE_PRIMERS}_${PRIMERS}_readsout.csv # columns: stage; sample_primers; read_group; primers; fwd_out; fwd_out

else 
    echo "PAIRED variable must be true or false"
    exit 1
fi 

