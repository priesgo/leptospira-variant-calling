# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs HaplotypeCaller for haploid organisms"
    echo "USAGE: haplotype_caller.sh INPUT_BAM OUTPUT_VCF REFERENCE"
    exit 1
fi

#Configuration
BASEDIR=$(dirname "$0")
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

# Haplotype caller variant calling pipeline
java -jar $GATK -T HaplotypeCaller -R $REFERENCE -I $INPUT_BAM --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $OUTPUT_VCF --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff -ploidy 1 -allowNonUniqueKmersInRef -nda -bamout $INPUT_BAM.hc_reassembly.bam -forceActive -disableOptimizations
