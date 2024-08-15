#!/bin/bash
set -e
set -u
## args are the following:
# $1 = reads_paths, one or two read paths separated by ;
# $2 = meta.for_primer_seq, aka. fwd primer sequence
# $3 = meta.rev_primer_seq, aka. rev primer sequence
# $4 = meta.pcr_primers, aka. name of PCR primer pair
# $5 = meta.target_gene, aka. name of target gene
# $6 = meta.sample_id, aka. sample ID
# $7 = meta.fcid, aka. flowcell ID
# $8 = params.primer_error_rate
# $9 = seq_type (illumina, nanopore, pacbio)
# $10 = paired (true, false)

### define variables with better names

READS_PATHS=$1
FOR_PRIMER_SEQ=$2
REV_PRIMER_SEQ=$3
PCR_PRIMERS=$4
TARGET_GENE=$5
SAMPLE_ID=$6
FCID=$7
PRIMER_ERROR_RATE=$8
SEQ_TYPE=$9
PAIRED=${10}


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
FWD_PRIMER=${FOR_PRIMER_SEQ/I/N}
REV_PRIMER=${REV_PRIMER_SEQ/I/N}

# reverse complement the primer sequences
FWD_PRIMER_RC=$(echo ${FWD_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)
REV_PRIMER_RC=$(echo ${REV_PRIMER} | \
    tr GATCRYMKSWHBVDNgatcrymkswhbvdn CTAGYRKMSWDVBHNctagyrkmswdvbhn | \
    rev)


### split reads, with conditional code based on single vs paired

if [ $PAIRED == "true" ]; then 

    ### paired-end trimming

    ## retain reads that contain primer sequence
    cutadapt \
        -e $PRIMER_ERROR_RATE \
        --action=retain \
        --match-read-wildcards \
        --report=minimal \
        --discard-untrimmed \
        --rename="{header}" \
        -g ^$FWD_PRIMER \
        -G ^$REV_PRIMER \
        -o ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R1.fastq.gz \
        -p ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R2.fastq.gz \
        $FWD_READS \
        $REV_READS

    ## count reads in input files
    # forward reads
    if [[ $FWD_READS == *.gz ]]; then
        R1_IN_LINES=$(zcat $FWD_READS | wc -l)
    else 
        R1_IN_LINES=$(cat $FWD_READS | wc -l)
    fi
    R1_IN=$(( $R1_IN_LINES / 4 ))

    # reverse reads
    if [[ $REV_READS == *.gz ]]; then
        R2_IN_LINES=$(zcat $REV_READS | wc -l)
    else 
        R2_IN_LINES=$(cat $REV_READS | wc -l)
    fi
    R2_IN=$(( $R2_IN_LINES / 4 ))

    # save as .csv 
    echo "input,$SAMPLE_ID,$FCID,$PCR_PRIMERS,$R1_IN,$R2_IN" > input_${SAMPLE_ID}_${PCR_PRIMERS}_readsin.csv # columns: stage; sample_id; fcid; pcr_primers; fwd_in; rev_in

    ## count reads in output files
    # forward reads
    if [ -f ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R1.fastq.gz ]; then
        R1_OUT_LINES=$(zcat ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R1.fastq.gz | wc -l)
        R1_OUT=$(( $R1_OUT_LINES / 4 ))
    else 
        R1_OUT=0
    fi

    # reverse reads
    if [ -f ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R2.fastq.gz ]; then
        R2_OUT_LINES=$(zcat ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R2.fastq.gz | wc -l)
        R2_OUT=$(( $R2_OUT_LINES / 4 ))
    else 
        R2_OUT=0
    fi

    # save as .csv 
    echo "split_loci,$SAMPLE_ID,$FCID,$PCR_PRIMERS,$R1_OUT,$R2_OUT" > split_loci_${SAMPLE_ID}_${PCR_PRIMERS}_readsout.csv # columns: stage; sample_id; fcid; pcr_primers; fwd_out; fwd_out

elif [ $PAIRED == "false" ]; then

    ### single-end trimming
    cutadapt \
        -e $PRIMER_ERROR_RATE \
        --action=retain \
        --match-read-wildcards \
        --report=minimal \
        --discard-untrimmed \
        --rename="{header}" \
        -a ${FWD_PRIMER}..${REV_PRIMER_RC} \
        --revcomp \
        -o ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R0.fastq.gz \
        $SINGLE_READS 

    ## count reads in input file
    # single reads
    if [[ $SINGLE_READS == *.gz ]]; then
        R0_IN_LINES=$(zcat $SINGLE_READS | wc -l)
    else 
        R0_IN_LINES=$(cat $SINGLE_READS | wc -l)
    fi
    R0_IN=$(( $R0_IN_LINES / 4 ))

    # save as .csv 
    echo "input,$SAMPLE_ID,$FCID,$PCR_PRIMERS,$R0_IN,$R0_IN" > input_${SAMPLE_ID}_${PCR_PRIMERS}_readsin.csv # columns: stage; sample_id; fcid; pcr_primers; fwd_in; rev_in

    ## count reads in output files
    # forward reads
    if [ -f ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R0.fastq.gz ]; then
        R0_OUT_LINES=$(zcat ${SAMPLE_ID}_${TARGET_GENE}_${PCR_PRIMERS}_R0.fastq.gz | wc -l)
        R0_OUT=$(( $R0_OUT_LINES / 4 ))
    else 
        R0_OUT=0
    fi

    # save as .csv 
    echo "split_loci,$SAMPLE_ID,$FCID,$PCR_PRIMERS,$R0_OUT,$R0_OUT" > split_loci_${SAMPLE_ID}_${PCR_PRIMERS}_readsout.csv # columns: stage; sample_id; fcid; pcr_primers; fwd_out; fwd_out

else 
    echo "PAIRED variable must be true or false"
    exit 1
fi 

