#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs Samtools variant calling"
    echo "USAGE: samtools_pileup.sh INPUT_BAM OUTPUT_VCF REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
PREFIX_LOCAL=`basename $1 .bam`
echo $INPUT_BAM
OUTPUT_VCF=$2
echo $OUTPUT_VCF
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
REFERENCE=$3
echo $REFERENCE

# Variant calling with samtools pileup
echo "Samtools variant calling"
$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE $INPUT_BAM | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_VCF
