# Parameters
SNPEFF=/opt/snpEff/snpEff.jar
SNPEFF_CONFIG=/opt/snpEff/snpEff.config

INPUT_VCF=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
SNPEFF_REFERENCE=$3
echo $SNPEFF_REFERENCE

# Translate NCBI identifiers into the sequence numbers used in SnpEff reference
# matching was confirmed using the number of bp of each chromosome
sed -i 's/NC_008510/1/g' $INPUT_VCF
sed -i 's/NC_008511/2/g' $INPUT_VCF
sed -i 's/NC_008508/1/g' $INPUT_VCF
sed -i 's/NC_008509/2/g' $INPUT_VCF

java -jar $SNPEFF eff $SNPEFF_REFERENCE $INPUT_VCF -c $SNPEFF_CONFIG -s $PREFIX_LOCAL.stats.html -canon > $OUTPUT_VCF
mv snpEff_genes.txt $PREFIX_LOCAL.genes.stats.txt

