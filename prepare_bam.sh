#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Preprocess BAM (clean, fix mate info, remove duplicates) and then realigns it around indels."
	echo "USAGE: prepare_bam.sh <INPUT_BAM> <OUTPUT_DIR> <REFERENCE>"
	exit 1
fi

# Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/config/config.sh

INPUT_BAM=$1
PREFIX=`basename $1`
# echo "$INPUT_BAM"
OUTPUT_DIR=$2
# echo "$OUTPUT_DIR"
REFERENCE=$3
# echo "$REFERENCE"
REFERENCE_BASENAME=`basename $REFERENCE`

echo "Preparing input BAM ..."
echo

# BAM preprocessing
echo "$BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_DIR/$PREFIX.preprocessed.bam"
echo
# source $BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_DIR/$PREFIX.preprocessed.bam

# Realignment around indels
echo "$BASEDIR/realignment/realign_bam.sh $OUTPUT_DIR/$PREFIX.preprocessed.bam $OUTPUT_DIR/$PREFIX.realigned.bam $REFERENCE"
echo
source $BASEDIR/realignment/realign_bam.sh $OUTPUT_DIR/$PREFIX.preprocessed.bam $OUTPUT_DIR/$PREFIX.realigned.bam $REFERENCE
echo
