# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs Samtools variant calling"
    echo "USAGE: samtools_pileup.sh INPUT_BAM OUTPUT_VCF REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_BAM
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE


# Variant calling with samtools pileup
$SAMTOOLS_HOME/samtools mpileup -6EOsuf $REFERENCE $INPUT_BAM | $BCFTOOLS_HOME/bcftools call -mv -Ov --ploidy 1 > $PREFIX_LOCAL.vcf
java -jar $GATK -T VariantAnnotator -I $INPUT_BAM -R $REFERENCE -V $PREFIX_LOCAL.vcf -o $OUTPUT_VCF --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
rm -f $PREFIX_LOCAL.vcf
rm -f $PREFIX_LOCAL.vcf
