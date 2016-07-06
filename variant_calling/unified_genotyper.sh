#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Runs UnifiedGenotyper for haploid organisms"
		echo "USAGE: unified_genotyper.sh INPUT_BAM OUTPUT_VCF REFERENCE"
    	exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/../config/config.sh

INPUT_BAM_LOCAL=$1
# echo $INPUT_BAM_LOCAL
OUTPUT_VCF_LOCAL=$2
# echo $OUTPUT_VCF_LOCAL
OUTPUT_DIR_LOCAL=$(dirname "$OUTPUT_VCF_LOCAL")
# echo $OUTPUT_VCF_LOCAL
# echo $OUTPUT_DIR_LOCAL
if [ ! -d "$OUTPUT_DIR_LOCAL" ]; then
	mkdir $OUTPUT_DIR_LOCAL
fi
REFERENCE_LOCAL=$3
# echo $REFERENCE_LOCAL

# Unified Genotyper variant calling pipeline
echo "Variant calling with GATK (UnifiedGenotyper) ..."
echo
echo "java -jar $GATK -T UnifiedGenotyper -R $REFERENCE_LOCAL -I $INPUT_BAM_LOCAL --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE -ploidy 1 -nda -o $OUTPUT_VCF_LOCAL --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T UnifiedGenotyper -R $REFERENCE_LOCAL -I $INPUT_BAM_LOCAL --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -baq RECALCULATE --downsampling_type NONE -ploidy 1 -nda -o $OUTPUT_VCF_LOCAL --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
echo
