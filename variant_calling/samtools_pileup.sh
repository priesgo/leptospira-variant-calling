#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Runs Samtools variant calling"
		echo "USAGE: samtools_pileup.sh INPUT_BAM OUTPUT_VCF REFERENCE"
		exit 1
fi

# Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

INPUT_BAM_ST=$1
PREFIX_ST=`basename $1 .bam`
# echo $INPUT_BAM_ST
# echo $PREFIX_ST
OUTPUT_VCF_ST=$2
OUTPUT_DIR_ST=$(dirname "$OUTPUT_VCF_ST")
# echo $OUTPUT_VCF_ST
# echo $OUTPUT_DIR_ST
if [ ! -d "$OUTPUT_DIR_ST" ]; then
  mkdir $OUTPUT_DIR_ST
fi
REFERENCE_ST=$3
# echo $REFERENCE_ST

# Variant calling with samtools pileup
echo "Variant calling with samtools ..."
echo
echo "$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE_ST $INPUT_BAM_ST | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR/$PREFIX_ST.tmp.vcf"
$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE_ST $INPUT_BAM_ST | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR_ST/$PREFIX_ST.tmp.vcf 
echo
echo "java -jar $GATK -T VariantAnnotator -I $INPUT_BAM_ST -R $REFERENCE_ST -V $OUTPUT_DIR_ST/$PREFIX_ST.tmp.vcf -o $OUTPUT_VCF_ST --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T VariantAnnotator -I $INPUT_BAM_ST -R $REFERENCE_ST -V $OUTPUT_DIR_ST/$PREFIX_ST.tmp.vcf -o $OUTPUT_VCF_ST --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff

rm -f $OUTPUT_DIR_ST/$PREFIX_ST.tmp.vcf
echo
