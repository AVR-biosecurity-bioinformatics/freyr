#!/bin/bash
set -e
set -u
## args are the following:
# $1 = primers
# $2 = read_group
# $3 = fasta (ASVs)
# $4 = blast_min_identity
# $5 = blast_min_coverage
# $6 = run_blast

### define variables with better names

PRIMERS=$1
READ_GROUP=$2
FASTA=$3
RUN_BLAST=$4


THREADS=$(nproc)

if [ $RUN_BLAST == "TRUE" ]; then

	# replace each space in sequence headers of fasta files with the string "!?!?"
	sed '/^>/ s/ /!?!?/g' $FASTA > query.fasta

	# run blast
	blastn \
		-query query.fasta \
		-db ${PRIMERS}.blast_db \
		-out blast.tsv \
		-outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs" \
		-num_threads $THREADS \
		-max_target_seqs 5

	rm -f query.fasta

else 
	
	touch blast.tsv

fi