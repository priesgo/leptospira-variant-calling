#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Preprocess BAM (clean, fix mate info, remove duplicates) and then realigns it around indels."
	echo "USAGE: prepare_bam.sh <INPUT_BAM> <OUTPUT_DIR> <REFERENCE>"
	exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/config/config.sh

INPUT_BAM_LOCAL=$1
PREFIX_LOCAL=`basename $1`
# echo "$INPUT_BAM"
OUTPUT_DIR_LOCAL=$2
# echo "$OUTPUT_DIR_LOCAL"
REFERENCE_LOCAL=$3
# echo "$REFERENCE"
# REFERENCE_BASENAME=`basename $REFERENCE_LOCAL`

echo "Preparing input BAM ..."
echo

# BAM preprocessingINPUT_BAM_LOCAL
echo "$BASEDIR_LOCAL/preprocessing/preprocess_bam.sh $INPUT_BAM_LOCAL $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.preprocessed.bam"
echo
source $BASEDIR_LOCAL/preprocessing/preprocess_bam.sh $INPUT_BAM_LOCAL $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.preprocessed.bam

# Realignment around indels
echo "$BASEDIR_LOCAL/realignment/realign_bam.sh $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.preprocessed.bam $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.realigned.bam $REFERENCE_LOCAL"
echo
source $BASEDIR_LOCAL/realignment/realign_bam.sh $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.preprocessed.bam $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.realigned.bam $REFERENCE_LOCAL
echo
