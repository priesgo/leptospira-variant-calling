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

INPUT_BAM=$1
PREFIX_LOCAL=`basename $1 .bam`
# echo $INPUT_BAM
# echo $PREFIX_LOCAL
OUTPUT_VCF=$2
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
# echo $OUTPUT_VCF
# echo $OUTPUT_DIR
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir $OUTPUT_DIR
fi
REFERENCE=$3
# echo $REFERENCE

# Variant calling with samtools pileup
echo "Variant calling with samtools ..."
echo
echo "$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE $INPUT_BAM | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf"
$SAMTOOLS_HOME/samtools mpileup --min-BQ 13 --adjust-MQ 50 --redo-BAQ --min-MQ 1 --illumina1.3+ --output-BP --output-MQ --uncompressed --fasta-ref $REFERENCE $INPUT_BAM | $BCFTOOLS_HOME/bcftools call --multiallelic-caller --variants-only --output-type v --ploidy 1 > $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf 
echo "java -jar $GATK -T VariantAnnotator -I $INPUT_BAM -R $REFERENCE -V $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf -o $OUTPUT_VCF --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff"
java -jar $GATK -T VariantAnnotator -I $INPUT_BAM -R $REFERENCE -V $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf -o $OUTPUT_VCF --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff

rm -f $OUTPUT_DIR/$PREFIX_LOCAL.tmp.vcf
echo
