#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Runs UnifiedGenotyper for haploid organisms"
		echo "USAGE: unified_genotyper.sh INPUT_BAM OUTPUT_VCF REFERENCE"
    	exit 1
fi

# Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

INPUT_BAM_UG=$1
# echo $INPUT_BAM_UG
OUTPUT_VCF_UG=$2
# echo $OUTPUT_VCF_UG
OUTPUT_DIR_UG=$(dirname "$OUTPUT_VCF_UG")
# echo $OUTPUT_VCF_UG
# echo $OUTPUT_DIR_UG
if [ ! -d "$OUTPUT_DIR_UG" ]; then
	mkdir $OUTPUT_DIR_UG
fi
REFERENCE_UG=$3
# echo $REFERENCE_UG

# Unified Genotyper variant calling pipeline
echo "Variant calling with GATK (UnifiedGenotyper) ..."
echo
echo "java -jar $GATK -T UnifiedGenotyper -R $REFERENCE_UG -I $INPUT_BAM_UG --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE -ploidy 1 -nda -o $OUTPUT_VCF_UG --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T UnifiedGenotyper -R $REFERENCE_UG -I $INPUT_BAM_UG --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE -ploidy 1 -nda -o $OUTPUT_VCF_UG --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
echo
