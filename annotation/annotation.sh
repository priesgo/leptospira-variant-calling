# Parameters
SNPEFF=/opt/snpEff/snpEff.jar
SNPEFF_CONFIG=/opt/snpEff/snpEff.config
GATK=/opt/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

INPUT_VCF=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_VCF
OUTPUT_VCF=$2
echo $OUTPUT_VCF
SNPEFF_REFERENCE=$3
echo $SNPEFF_REFERENCE
REFERENCE=$4
echo $REFERENCE

REFERENCE_BASENAME=`basename $REFERENCE`

if ["$REFERENCE_BASENAME"=="Lb.Hardjo.L550.fasta"] 
then
	CHR1="NC_008508";
        CHR2="NC_008509";
else #if ["$REFERENCE_BASENAME" == "Lb.Hardjo.JB197.fasta"] then
	CHR1="NC_008510";
        CHR2="NC_008511";
#else
#	echo "Reference not recognised as Leptospira";
#	exit 1;
fi


# Translate NCBI identifiers into the sequence numbers used in SnpEff reference
# matching was confirmed using the number of bp of each chromosome
cp $INPUT_VCF $PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008510/1/g' $PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008511/2/g' $PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008508/1/g' $PREFIX_LOCAL.snpeff_ref.vcf
sed -i 's/NC_008509/2/g' $PREFIX_LOCAL.snpeff_ref.vcf

java -jar $SNPEFF eff $SNPEFF_REFERENCE $PREFIX_LOCAL.snpeff_ref.vcf -c $SNPEFF_CONFIG -s $PREFIX_LOCAL.stats.html -canon -o gatk -hgvs > $PREFIX_LOCAL.annotated_all_transcripts.vcf
mv snpEff_genes.txt $PREFIX_LOCAL.genes.stats.txt

# Pastes the VCF header
cat $PREFIX_LOCAL.annotated_all_transcripts.vcf | grep '#' > $PREFIX_LOCAL.ncbi_ref.vcf
cat $PREFIX_LOCAL.annotated_all_transcripts.vcf | grep -v '#' | awk -v CHR1="$CHR1" -v CHR2="$CHR2" 'BEGIN {FS = "\t";OFS = "\t"}{if ($1=="1") {$1=CHR1; print;} else {$1=CHR2; print;}}' >> $PREFIX_LOCAL.ncbi_ref.vcf

java -jar $GATK -T VariantAnnotator -R $REFERENCE -A SnpEff --variant $INPUT_VCF --snpEffFile $PREFIX_LOCAL.ncbi_ref.vcf -L $INPUT_VCF -o $OUTPUT_VCF

