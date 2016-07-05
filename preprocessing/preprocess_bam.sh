#!/bin/bash

# Check input parameters
if [ $# -ne 2 ]
	then
		echo "Runs Picard preprocessing on a BAM file"
		echo "USAGE: preprocess_bam.sh INPUT_BAM OUTPUT_BAM"
    	exit 1
fi

# Configuration
SCRIPT_LOCAL_2=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL_2=$(dirname "$SCRIPT_LOCAL_2")
source $BASEDIR_LOCAL_2/../config/config.sh

INPUT_BAM_LOCAL_2=$1
PREFIX_LOCAL_2=`basename $1 .bam`
OUTPUT_BAM_LOCAL_2=$2
OUTPUT_DIR_LOCAL_2=$(dirname "$OUTPUT_BAM_LOCAL_2")

# Preprocessing pipeline
echo "Picard CleanSam"
echo
java -jar $PICARD CleanSam I=$INPUT_BAM_LOCAL_2 O=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.cleaned.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
echo
echo "Picard AddOrReplaceReadGroups"
echo
java -jar $PICARD AddOrReplaceReadGroups I=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.cleaned.bam O=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.bam LB=Library PL=Illumina PU=Barcode SM=$PREFIX_LOCAL_2  VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
echo
echo "Picard SortSam"
echo
java -jar $PICARD SortSam I=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.bam O=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.sorted.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=$TMP_DIR
echo
echo "Picard MarkDuplicates"
echo
java -jar $PICARD MarkDuplicates I=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.sorted.bam O=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.dedupped.bam METRICS_FILE=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.dedupped.bam.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR PROGRAM_RECORD_ID=null
echo
echo "Picard FixMateInformation"
echo
java -jar $PICARD FixMateInformation I=$OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.dedupped.bam O=$OUTPUT_BAM_LOCAL_2 VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=false TMP_DIR=$TMP_DIR
echo

rm -f $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.cleaned.ba*
rm -f $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.ba*
rm -f $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.readgroups.sorted.ba*
rm -f $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.dedupped.ba*
