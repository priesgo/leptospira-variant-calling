#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Recalibrates mapping qualities based on covariates"
    echo "USAGE: recalibrate_mapping_qualities.sh INPUT_BAM OUTPUT_BAM REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
PREFIX_LOCAL=`basename $1 .bam`
INPUT_BAM=$1
echo $INPUT_BAM
OUTPUT_BAM=$2
echo $OUTPUT_BAM
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
REFERENCE=$3
echo $REFERENCE

# Mapping quality recalibration pipeline
#java -Xmx1g -jar $GATK -T BaseRecalibrator -I $INPUT_BAM -R $REFERENCE -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate  --out $PREFIX_LOCAL.recal_data.csv --filter_mismatching_base_and_quals --allow_potentially_misencoded_quality_scores -solid_nocall_strategy LEAVE_READ_UNRECALIBRATED --disable_indel_quals

java -Xmx1g -jar $GATK -T PrintReads -I $INPUT_BAM -R $REFERENCE -BQSR $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp --out $OUTPUT_BAM --filter_mismatching_base_and_quals --allow_potentially_misencoded_quality_scores

java -jar $GATK -T AnalyzeCovariates -R $REFERENCE -BQSR $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp -plots $OUTPUT_DIR/$PREFIX_LOCAL.BQSR.pdf


rm -f $OUTPUT_DIR/$PREFIX_LOCAL.recal_data.grp
