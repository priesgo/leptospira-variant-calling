# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Filters a VCF file using recommended hard thresholds"
    echo "USAGE: variant_filtering.sh INPUT_VCF OUTPUT_VCF REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Parameters
INPUT_VCF=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE

# Variant filtering pipeline for SNPs
java -jar $GATK -T SelectVariants -R $REFERENCE -V $INPUT_VCF -selectType SNP -o $PREFIX_LOCAL.raw_snps.vcf
java -jar $GATK -T VariantFiltration -R $REFERENCE -V $PREFIX_LOCAL.raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filterName "FP" -o $PREFIX_LOCAL.filtered_snps.vcf

# Variant filtering pipeline for indels
java -jar $GATK -T SelectVariants -R $REFERENCE -V $INPUT_VCF -selectType INDEL -o $PREFIX_LOCAL.raw_indels.vcf
java -jar $GATK -T VariantFiltration -R $REFERENCE -V $PREFIX_LOCAL.raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0" --filterName "FP" -o $PREFIX_LOCAL.filtered_indels.vcf

# Combine all variants 
# --filteredAreUncalled
java -jar $GATK -T CombineVariants -R $REFERENCE -o $OUTPUT_VCF -V $PREFIX_LOCAL.filtered_snps.vcf -V $PREFIX_LOCAL.filtered_indels.vcf --genotypemergeoption UNSORTED

rm -f $PREFIX_LOCAL.filtered_snps.vcf
rm -f $PREFIX_LOCAL.filtered_snps.vcf.idx
rm -f $PREFIX_LOCAL.filtered_indels.vcf
rm -f $PREFIX_LOCAL.filtered_indels.vcf.idx
rm -f $PREFIX_LOCAL.raw_indels.vcf
rm -f $PREFIX_LOCAL.raw_indels.vcf.idx
rm -f $PREFIX_LOCAL.raw_snps.vcf
rm -f $PREFIX_LOCAL.raw_snps.vcf.idx

