#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs HaplotypeCaller for haploid organisms"
    echo "USAGE: haplotype_caller.sh INPUT_BAM OUTPUT_VCF REFERENCE"
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
PREFIX_LOCAL=`basename $1 .bam`
OUTPUT_VCF=$2
echo $OUTPUT_VCF
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
REFERENCE=$3
echo $REFERENCE

# Runs BAQ with PrintReads as HaplotypeCaller does not support it on the fly as UnifiedGenotyper does
echo "GATK PrintReads to calculate BAQ"
java -jar $GATK -T PrintReads -R $REFERENCE -I $INPUT_BAM -baq RECALCULATE -o $OUTPUT_DIR/$PREFIX_LOCAL.baq.bam

# Haplotype caller variant calling pipeline
echo "GATK HaplotypeCaller"
java -jar $GATK -T HaplotypeCaller -R $REFERENCE -I $OUTPUT_DIR/$PREFIX_LOCAL.baq.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 -ploidy 1 -nda -allowNonUniqueKmersInRef -bamout $INPUT_BAM.hc_reassembly.bam -o $OUTPUT_VCF --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
# Use this parameters to create all haplotypes -forceActive -disableOptimizations

rm -f $OUTPUT_DIR/$PREFIX_LOCAL.baq.ba*
