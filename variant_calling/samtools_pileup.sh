#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Runs Samtools variant calling"
		echo "USAGE: samtools_pileup.sh INPUT_BAM OUTPUT_VCF REFERENCE"
		exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/../config/config.sh

INPUT_BAM_LOCAL=$1
PREFIX_LOCAL=`basename $1 .bam`
# echo $INPUT_BAM_LOCAL
# echo $PREFIX_LOCAL
OUTPUT_VCF_LOCAL=$2
OUTPUT_DIR_LOCAL=$(dirname "$OUTPUT_VCF_LOCAL")
# echo $OUTPUT_VCF_LOCAL
# echo $OUTPUT_DIR_LOCAL
if [ ! -d "$OUTPUT_DIR_LOCAL" ]; then
  mkdir $OUTPUT_DIR_LOCAL
fi
REFERENCE_LOCAL=$3
# echo $REFERENCE_LOCAL

# Variant calling with samtools pileup
echo "Variant calling with samtools ..."
echo
echo "$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE_LOCAL $INPUT_BAM_LOCAL | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf"
$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE_LOCAL $INPUT_BAM_LOCAL | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.tmp.vcf 
echo
echo "java -jar $GATK -T VariantAnnotator -I $INPUT_BAM_LOCAL -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.tmp.vcf -o $OUTPUT_VCF_LOCAL --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T VariantAnnotator -I $INPUT_BAM_LOCAL -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.tmp.vcf -o $OUTPUT_VCF_LOCAL --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
echo
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.tmp.vcf
