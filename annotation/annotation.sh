#!/bin/bash
# Check input parameters
if [ $# -ne 4 ]
  then
	echo "Annotates VCF file for Leptospira borgspetersenii L550 and JB197 strains"
    echo "USAGE: annotation.sh INPUT_VCF OUTPUT_VCF SNPEFF_REFERENCE REFERENCE"
    echo "Alternative for SNPEFF_REFERENCE are: Leptospira_borgpetersenii_serovar_Hardjo_bovis_L550_uid58507 and Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

INPUT_VCF=$1
PREFIX_LOCAL=`basename $1 .vcf`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
SNPEFF_REFERENCE=$3
echo $SNPEFF_REFERENCE
REFERENCE=$4
echo $REFERENCE

REFERENCE_BASENAME=`basename $REFERENCE`
L550="Lb.Hardjo.L550.fasta"

if [ "$REFERENCE_BASENAME" == "$L550" ]
then
	CHR1="NC_008508";
        CHR2="NC_008509";
else
	CHR1="NC_008510";
        CHR2="NC_008511";
fi


# Translate NCBI identifiers into the sequence numbers used in SnpEff reference
# matching was confirmed using the number of bp of each chromosome
cp $INPUT_VCF $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008510/1/g' $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008511/2/g' $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008508/1/g' $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008509/2/g' $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf

# -o gatk required to select the highest impact transcript using olf format
# java -jar $SNPEFF eff $SNPEFF_REFERENCE $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf -c $SNPEFF_CONFIG -s $OUTPUT_DIR/$PREFIX_LOCAL.stats.html -canon -hgvs -ud 0 -onlyProtein > $OUTPUT_DIR/$PREFIX_LOCAL.annotated_all_transcripts.vcf

java -jar /opt/snpEff/snpEff.jar eff $SNPEFF_REFERENCE $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf -c $SNPEFF_CONFIG -s $OUTPUT_DIR/$PREFIX_LOCAL.stats.html -canon -hgvs -ud 0 -onlyProtein > $OUTPUT_DIR/$PREFIX_LOCAL.annotated_all_transcripts.vcf
# mv $OUTPUT_DIR/snpEff_genes.txt $OUTPUT_DIR/$PREFIX_LOCAL.genes.stats.txt

# Pastes the VCF header
cat $OUTPUT_DIR/$PREFIX_LOCAL.annotated_all_transcripts.vcf | grep '#' > $OUTPUT_VCF
cat $OUTPUT_DIR/$PREFIX_LOCAL.annotated_all_transcripts.vcf | grep -v '#' | awk -v CHR1="$CHR1" -v CHR2="$CHR2" 'BEGIN {FS = "\t";OFS = "\t"}{if ($1=="1") {$1=CHR1; print;} else {$1=CHR2; print;}}' >> $OUTPUT_VCF

rm -f $OUTPUT_DIR/$PREFIX_LOCAL.snpeff_ref.vcf*
rm -f $OUTPUT_DIR/$PREFIX_LOCAL.annotated_all_transcripts.vcf
