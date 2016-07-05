#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Filters a VCF file using recommended hard thresholds"
		echo "USAGE: variant_filtering.sh INPUT_VCF OUTPUT_VCF REFERENCE"
		exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/../config/config.sh

# Parameters
INPUT_VCF_LOCAL=$1
PREFIX_LOCAL=`basename $1 .vcf`
# echo $INPUT_VCF_LOCAL
OUTPUT_VCF_LOCAL=$2
# echo $OUTPUT_VCF_LOCAL
OUTPUT_DIR_LOCAL=$(dirname "$OUTPUT_VCF_LOCAL")
REFERENCE_LOCAL=$3
# echo $REFERENCE_LOCAL

echo "Filtering GATK variants ..."
echo

# Variant filtering pipeline for SNPs
echo "java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V $INPUT_VCF_LOCAL -selectType SNP -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf"
java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V $INPUT_VCF_LOCAL -selectType SNP -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf
echo
echo "java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf --filterExpression \"QD < 2.0\" --filterName \"QD\" --filterExpression \"SOR > 6.0\" --filterName \"SOR\" --filterExpression \"QUAL < 5\" --filterName \"QUAL\" --filterExpression \"DP < 3\" --filterName \"DP\" -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf"
java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "SOR > 6.0" --filterName "SOR" --filterExpression "QUAL < 5" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf
echo

# Variant filtering pipeline for indels
echo "java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V $INPUT_VCF_LOCAL -selectType INDEL -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf"
java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V $INPUT_VCF_LOCAL -selectType INDEL -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf
echo
echo "java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf --filterExpression \"QD < 2.0\" --filterName \"QD\" --filterExpression \"SOR > 10.0\" --filterName \"SOR\" --filterExpression \"QUAL < 5\" --filterName \"QUAL\" --filterExpression \"DP < 3\" --filterName \"DP\" -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf"
java -jar $GATK -T VariantFiltration -R $REFERENCE_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "SOR > 10.0" --filterName "SOR" --filterExpression "QUAL < 5" --filterName "QUAL" --filterExpression "DP < 3" --filterName "DP" -o $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf
echo

# Combine all variants 
# --filteredAreUncalled
echo "java -jar $GATK -T CombineVariants -R $REFERENCE_LOCAL -o $OUTPUT_VCF_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf --genotypemergeoption UNSORTED"
java -jar $GATK -T CombineVariants -R $REFERENCE_LOCAL -o $OUTPUT_VCF_LOCAL -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf -V $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf --genotypemergeoption UNSORTED
echo

rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_snps.vcf.idx
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.filtered_indels.vcf.idx
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_indels.vcf.idx
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf
rm -f $OUTPUT_DIR_LOCAL/$PREFIX_LOCAL.raw_snps.vcf.idx

