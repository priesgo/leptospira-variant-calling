#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Preprocess BAM (clean, fix mate info, remove duplicates) and then realigns it around indels."
	echo "USAGE: prepare_bam.sh <INPUT_BAM> <OUTPUT_FOLDER> <REFERENCE>"
	exit 1
fi

#mkdir logs
#> logs/$PREFIX.log
#2> logs/$PREFIX.log

# Parameters
SCRIPT=$(readlink -f "$BASH_SOURCE")
PREPARE_BAM_BASEDIR=$(dirname "$SCRIPT")
INPUT_BAM=$1
PREFIX=`basename $1`
# echo "$INPUT_BAM"
OUTPUT_FOLDER=$2
# echo "$OUTPUT_FOLDER"
REFERENCE=$3
# echo "$REFERENCE"
REFERENCE_BASENAME=`basename $REFERENCE`

echo "Preparing input BAM ..."
echo

# BAM preprocessing
echo "$PREPARE_BAM_BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_FOLDER/$PREFIX.preprocessed.bam"
source $PREPARE_BAM_BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_FOLDER/$PREFIX.preprocessed.bam

# Realignment around indels
echo "$PREPARE_BAM_BASEDIR/realignment/realign_bam.sh $OUTPUT_FOLDER/$PREFIX.preprocessed.bam $OUTPUT_FOLDER/$PREFIX.realigned.bam $REFERENCE"
source $PREPARE_BAM_BASEDIR/realignment/realign_bam.sh $OUTPUT_FOLDER/$PREFIX.preprocessed.bam $OUTPUT_FOLDER/$PREFIX.realigned.bam $REFERENCE
echo
