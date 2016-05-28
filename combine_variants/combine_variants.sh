# Check input parameters
if [ $# -ne 6 ]
  then
	echo "Combines the variant calls of HaplotypeCaller, UnifiedGenotyper and samtools in the intersection and union sets"
    echo "USAGE: combine_variants.sh HC_INPUT_VCF UG_INPUT_VCF ST_INPUT_VCF UNION_OUTPUT_VCF INTERSECTION_OUTPUT_VCF REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

HC_INPUT_VCF=$1
#PREFIX_LOCAL=`basename $1`
echo $HC_INPUT_VCF
UG_INPUT_VCF=$2
echo $UG_INPUT_VCF
ST_INPUT_VCF=$3
echo $ST_INPUT_VCF
UNION_VCF=$4
echo $UNION_VCF
INTERSECT_VCF=$5
echo $INTERSECT_VCF
REFERENCE=$6
echo $REFERENCE

# Variant filtering pipeline for SNPs
java -jar $GATK -T CombineVariants -R $REFERENCE -V:HC $HC_INPUT_VCF -V:UG $UG_INPUT_VCF -V:ST $ST_INPUT_VCF -o $UNION_VCF --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED --genotypemergeoption UNIQUIFY
java -jar $GATK -T SelectVariants -R $REFERENCE -V:variant $UNION_VCF -select 'set == "Intersection";' -o $INTERSECT_VCF --setFilteredGtToNocall

# Calculates statistics
java -jar $GATK -T VariantEval -R $REFERENCE --eval $UNION_VCF -o $UNION_VCF.gatkreport
java -jar $GATK -T VariantEval -R $REFERENCE --eval $INTERSECT_VCF -o $INTERSECT_VCF.gatkreport
