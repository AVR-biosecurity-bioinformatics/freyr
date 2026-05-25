#!/bin/bash
set -e
set -u
## args are the following:
# $1 = primers
# $2 = ref_fasta

### define variables with better names

PRIMERS=$1
REF_FASTA=$2

# check if ref_fasta is gzip compressed
if [[ $REF_FASTA == *.gz ]]; then
    gunzip -c $REF_FASTA > ref.fasta
else 
    cat $REF_FASTA > ref.fasta
fi

# replace each space in sequence headers of fasta file with the string "!?!?"
sed '/^>/ s/ /!?!?/g' ref.fasta > ref_converted.fasta

# create database
makeblastdb \
	-in ref_converted.fasta \
	-input_type fasta \
	-out ${PRIMERS}.blast_db \
	-dbtype nucl

# remove unneeded files
rm -f ref.fasta
rm -f ref_converted.fasta