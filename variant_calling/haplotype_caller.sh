#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Runs HaplotypeCaller for haploid organisms"
		echo "USAGE: haplotype_caller.sh INPUT_BAM OUTPUT_VCF REFERENCE"
		exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/../config/config.sh

INPUT_BAM_LOCAL=$1
PREFIX_LOCAL=`basename $1 .bam`
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

# Runs BAQ with PrintReads as HaplotypeCaller does not support it on the fly as UnifiedGenotyper does
echo "Variant calling with GATK (HaplotypeCaller) ..."
echo
echo "Running PrintReads to calculate BAQ"
echo
echo "java -jar $GATK -T PrintReads -R $REFERENCE_LOCAL -I $INPUT_BAM_LOCAL -baq RECALCULATE -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.baq.bam"
java -jar $GATK -T PrintReads -R $REFERENCE_LOCAL -I $INPUT_BAM_LOCAL -baq RECALCULATE -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.baq.bam
echo
# Haplotype caller variant calling pipeline
echo "Running HaplotypeCaller"
echo
echo "java -jar $GATK -T HaplotypeCaller -R $REFERENCE_LOCAL -I $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.baq.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 --downsampling_type NONE -ploidy 1 -nda -allowNonUniqueKmersInRef -bamout $INPUT_BAM_LOCAL.hc_reassembly.bam -o $OUTPUT_VCF_LOCAL --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T HaplotypeCaller -R $REFERENCE_LOCAL -I $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.baq.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 --downsampling_type NONE -ploidy 1 -nda -allowNonUniqueKmersInRef -bamout $INPUT_BAM_LOCAL.hc_reassembly.bam -o $OUTPUT_VCF_LOCAL --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
# Use this parameters to create all haplotypes -forceActive -disableOptimizations

# rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.baq.ba*
echo
