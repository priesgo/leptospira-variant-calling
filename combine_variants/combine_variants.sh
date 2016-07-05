#!/bin/bash
# Check input parameters
if [ $# -ne 6 ]
	then
		echo "Combines the variant calls of HaplotypeCaller, UnifiedGenotyper and samtools in the intersection and union sets"
		echo "USAGE: combine_variants.sh <HC_INPUT_VCF> <UG_INPUT_VCF> <ST_INPUT_VCF> <UNION_OUTPUT_VCF> <INTERSECTION_OUTPUT_VCF> <REFERENCE>"
		exit 1
fi

# Configuration
SCRIPT_LOCAL=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL=$(dirname "$SCRIPT_LOCAL")
source $BASEDIR_LOCAL/../config/config.sh

HC_INPUT_VCF=$1
# echo $HC_INPUT_VCF
UG_INPUT_VCF=$2
# echo $UG_INPUT_VCF
ST_INPUT_VCF=$3
# echo $ST_INPUT_VCF
UNION_VCF=$4
# echo $UNION_VCF
INTERSECT_VCF=$5
# echo $INTERSECT_VCF
REFERENCE_LOCAL=$6
# echo $REFERENCE_LOCAL

# Combine the variants and ouput the union and intersection files
echo "Combining variants ..."
echo
echo "java -jar $GATK -T CombineVariants -R $REFERENCE_LOCAL -V:HC $HC_INPUT_VCF -V:UG $UG_INPUT_VCF -V:ST $ST_INPUT_VCF -o $UNION_VCF --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED --genotypemergeoption UNIQUIFY"
java -jar $GATK -T CombineVariants -R $REFERENCE_LOCAL -V:HC $HC_INPUT_VCF -V:UG $UG_INPUT_VCF -V:ST $ST_INPUT_VCF -o $UNION_VCF --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED --genotypemergeoption UNIQUIFY
echo
echo "java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V:variant $UNION_VCF -select \'set == \"Intersection\";\' -o $INTERSECT_VCF --setFilteredGtToNocall"
java -jar $GATK -T SelectVariants -R $REFERENCE_LOCAL -V:variant $UNION_VCF -select 'set == "Intersection";' -o $INTERSECT_VCF --setFilteredGtToNocall
echo

# Calculates statistics
echo "Calculating stats and writing reports ..."
echo
echo "java -jar $GATK -T VariantEval -R $REFERENCE_LOCAL --eval $UNION_VCF -o $UNION_VCF.gatkreport"
java -jar $GATK -T VariantEval -R $REFERENCE_LOCAL --eval $UNION_VCF -o $UNION_VCF.gatkreport
echo
echo "java -jar $GATK -T VariantEval -R $REFERENCE_LOCAL --eval $INTERSECT_VCF -o $INTERSECT_VCF.gatkreport"
java -jar $GATK -T VariantEval -R $REFERENCE_LOCAL --eval $INTERSECT_VCF -o $INTERSECT_VCF.gatkreport
echo
