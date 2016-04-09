# Parameters
GATK=/opt/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
PREFIX=`basename $1`
echo $INPUT_BAM
OUTPUT_BAM=$2
echo $OUTPUT_BAM
REFERENCE=$3
echo $REFERENCE

# Preprocessing pipeline
java -Xmx1g -jar $GATK -T RealignerTargetCreator -R $REFERENCE -I $INPUT_BAM -o $PREFIX.intervals --filter_mismatching_base_and_quals 
java -Xmx4g -Djava.io.tmpdir=. -jar $GATK -T IndelRealigner -I $INPUT_BAM -R $REFERENCE -targetIntervals $PREFIX.intervals -o $OUTPUT_BAM --filter_mismatching_base_and_quals --consensusDeterminationModel USE_READS

rm -f $PREFIX.intervals
