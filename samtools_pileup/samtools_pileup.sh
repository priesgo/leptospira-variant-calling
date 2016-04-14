# Parameters
GATK=/opt/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
PREFIX_LOCAL=`basename $1`
echo $INPUT_BAM
OUTPUT_VCF=$2
echo $OUTPUT_VCF
REFERENCE=$3
echo $REFERENCE


# Variant calling with samtools pileup
samtools mpileup -DSugBf $REFERENCE $INPUT_BAM > $PREFIX_LOCAL.tmp.bcf
bcftools view -g $PREFIX_LOCAL.tmp.bcf > $PREFIX_LOCAL.bcf
bcftools view $PREFIX_LOCAL.bcf | vcfutils.pl varFilter - > $OUTPUT_VCF

rm -f $PREFIX_LOCAL.tmp.bcf
rm -f $PREFIX_LOCAL.bcf
