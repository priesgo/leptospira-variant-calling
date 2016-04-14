# Parameters
GATK=/opt/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

INPUT_VCF=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE

# Variant filtering pipeline for SNPs
java -jar $GATK -T SelectVariants -R $REFERENCE -V $INPUT_VCF -selectType SNP -o $PREFIX_LOCAL.raw_snps.vcf
java -jar $GATK -T VariantFiltration -R $REFERENCE -V $PREFIX_LOCAL.raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "FP" -o $PREFIX_LOCAL.filtered_snps.vcf

# Variant filtering pipeline for indels
java -jar $GATK -T SelectVariants -R $REFERENCE -V $INPUT_VCF -selectType INDEL -o $PREFIX_LOCAL.raw_indels.vcf
java -jar $GATK -T VariantFiltration -R $REFERENCE -V $PREFIX_LOCAL.raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "FP" -o $PREFIX_LOCAL.filtered_indels.vcf

# Combine all variants
java -jar $GATK -T CombineVariants -R $REFERENCE -o $OUTPUT_VCF -V $PREFIX_LOCAL.filtered_snps.vcf -V $PREFIX_LOCAL.filtered_indels.vcf -genotypeMergeOptions UNIQUIFY

rm -f $PREFIX_LOCAL.filtered_snps.vcf
rm -f $PREFIX_LOCAL.filtered_snps.vcf.idx
rm -f $PREFIX_LOCAL.filtered_indels.vcf
rm -f $PREFIX_LOCAL.filtered_indels.vcf.idx
rm -f $PREFIX_LOCAL.raw_indels.vcf
rm -f $PREFIX_LOCAL.raw_indels.vcf.idx
rm -f $PREFIX_LOCAL.raw_snps.vcf
rm -f $PREFIX_LOCAL.raw_snps.vcf.idx

