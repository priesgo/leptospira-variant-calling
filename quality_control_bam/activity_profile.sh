#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs quality controls on BAM file"
    echo "USAGE: quality_control_bam.sh <INPUT_BAM> <REFERENCE> <OUTPUT_FOLDER>"
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
REFERENCE=$2
echo $REFERENCE
OUTPUT_FOLDER=$3
echo $OUTPUT_FOLDER


echo "Run activity profile"
java -Xmx1g -jar $GATK -T FindCoveredIntervals --input_file $INPUT_BAM --reference_sequence $REFERENCE --activityProfileOut $OUTPUT_FOLDER/$PREFIX_LOCAL.activityprofile