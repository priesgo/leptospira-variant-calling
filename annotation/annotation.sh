# Parameters
SNPEFF=/opt/snpEff/snpEff.jar
SNPEFF_CONFIG=/opt/snpEff/snpEff.config

INPUT_VCF=$1
PREFIX=`basename $1`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE

java -jar $SNPEFF eff $REFERENCE $INPUT_VCF -c $SNPEFF_CONFIG > $OUTPUT_VCF

