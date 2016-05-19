# Parameters
source ../config/config.sh

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
PREFIX_LOCAL=`basename $1`
INPUT_BAM=$1
echo $INPUT_BAM
OUTPUT_BAM=$2
echo $OUTPUT_BAM
REFERENCE=$3
echo $REFERENCE

# Mapping quality recalibration pipeline
#java -Xmx1g -jar $GATK -T BaseRecalibrator -I $INPUT_BAM -R $REFERENCE -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate  --out $PREFIX_LOCAL.recal_data.csv --filter_mismatching_base_and_quals --allow_potentially_misencoded_quality_scores -solid_nocall_strategy LEAVE_READ_UNRECALIBRATED --disable_indel_quals

java -Xmx1g -jar $GATK -T PrintReads -I $INPUT_BAM -R $REFERENCE -BQSR $PREFIX_LOCAL.recal_data.grp --out $OUTPUT_BAM --filter_mismatching_base_and_quals --allow_potentially_misencoded_quality_scores

java -jar $GATK -T AnalyzeCovariates -R $REFERENCE -BQSR $PREFIX_LOCAL.recal_data.grp -plots $PREFIX_LOCAL.BQSR.pdf


rm -f $PREFIX_LOCAL.recal_data.grp
