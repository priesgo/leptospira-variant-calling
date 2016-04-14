# Parameters
DATA_FOLDER=/data/Leptospira/mapping/
PICARD=/opt/picard-tools-2.2.1/picard.jar
TMP_DIR=.

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_BAM
OUTPUT_BAM=$2
echo $OUTPUT_BAM

# Preprocessing pipeline
java -jar $PICARD CleanSam I=$INPUT_BAM O=$PREFIX_LOCAL.cleaned.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
java -jar $PICARD AddOrReplaceReadGroups I=$PREFIX_LOCAL.cleaned.bam O=$PREFIX_LOCAL.readgroups.bam LB=Library PL=Illumina PU=Barcode SM=$PREFIX_LOCAL  VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR
java -jar $PICARD MarkDuplicates I=$PREFIX_LOCAL.readgroups.bam O=$PREFIX_LOCAL.dedupped.bam METRICS_FILE=$PREFIX_LOCAL.dedupped.bam.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TMP_DIR PROGRAM_RECORD_ID=null
java -jar $PICARD FixMateInformation I=$PREFIX_LOCAL.dedupped.bam O=$PREFIX_LOCAL.fixed.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=false TMP_DIR=$TMP_DIR
java -jar $PICARD SortSam I=$PREFIX_LOCAL.fixed.bam O=$OUTPUT_BAM VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=$TMP_DIR
rm -f $PREFIX_LOCAL.cleaned.ba*
rm -f $PREFIX_LOCAL.readgroups.ba*
rm -f $PREFIX_LOCAL.dedupped.ba*
rm -f $PREFIX_LOCAL.fixed.ba*
