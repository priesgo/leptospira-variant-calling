#!/bin/bash
# Check input parameters
if [ $# -ne 2 ]
  then
	echo "Runs Picard preprocessing on a BAM file"
    echo "USAGE: preprocess_bam.sh INPUT_BAM OUTPUT_BAM"
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
OUTPUT_BAM=$2
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
echo $OUTPUT_BAM

# Preprocessing pipeline
java -jar $PICARD CleanSam I=$INPUT_BAM O=$OUTPUT_DIR/$PREFIX_LOCAL.cleaned.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
java -jar $PICARD AddOrReplaceReadGroups I=$OUTPUT_DIR/$PREFIX_LOCAL.cleaned.bam O=$OUTPUT_DIR/$PREFIX_LOCAL.readgroups.bam LB=Library PL=Illumina PU=Barcode SM=$PREFIX_LOCAL  VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
java -jar $PICARD MarkDuplicates I=$OUTPUT_DIR/$PREFIX_LOCAL.readgroups.bam O=$OUTPUT_DIR/$PREFIX_LOCAL.dedupped.bam METRICS_FILE=$OUTPUT_DIR/$PREFIX_LOCAL.dedupped.bam.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR PROGRAM_RECORD_ID=null
java -jar $PICARD FixMateInformation I=$OUTPUT_DIR/$PREFIX_LOCAL.dedupped.bam O=$OUTPUT_DIR/$PREFIX_LOCAL.fixed.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=false TMP_DIR=$TMP_DIR
java -jar $PICARD SortSam I=$OUTPUT_DIR/$PREFIX_LOCAL.fixed.bam O=$OUTPUT_BAM VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=$TMP_DIR
rm -f $OUTPUT_DIR/$PREFIX_LOCAL.cleaned.ba*
rm -f $OUTPUT_DIR/$PREFIX_LOCAL.readgroups.ba*
rm -f $OUTPUT_DIR/$PREFIX_LOCAL.dedupped.ba*
rm -f $OUTPUT_DIR/$PREFIX_LOCAL.fixed.ba*
