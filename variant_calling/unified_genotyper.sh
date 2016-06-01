#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs UnifiedGenotyper for haploid organisms"
    echo "USAGE: unified_genotyper.sh INPUT_BAM OUTPUT_VCF REFERENCE"
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
echo $INPUT_BAM
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE

# Haplotype caller variant calling pipeline
"GATK UnifiedGenotyper"
java -jar $GATK -T UnifiedGenotyper -R $REFERENCE -I $INPUT_BAM --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE -ploidy 1 -nda -o $OUTPUT_VCF --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff