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
SAMTOOLS_HOME=/home/priesgo/samtools-0.1.19
BCFTOOLS_HOME=$SAMTOOLS_HOME/bcftools
VCFUTILS_HOME=/home/priesgo/bcbio/bin


# Variant calling with samtools pileup
$SAMTOOLS_HOME/samtools mpileup -DSugBf $REFERENCE $INPUT_BAM > $PREFIX_LOCAL.tmp.bcf
$BCFTOOLS_HOME/bcftools view -bvcg $PREFIX_LOCAL.tmp.bcf > $PREFIX_LOCAL.bcf
$BCFTOOLS_HOME/bcftools index $PREFIX_LOCAL.bcf
$BCFTOOLS_HOME/bcftools view $PREFIX_LOCAL.bcf | $VCFUTILS_HOME/vcfutils.pl varFilter - > $PREFIX_LOCAL.vcf

rm -f $PREFIX_LOCAL.tmp.bcf
rm -f $PREFIX_LOCAL.bcf

java -jar $GATK -T VariantAnnotator -R $REFERENCE -V $PREFIX_LOCAL.vcf -o $OUTPUT_VCF --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation MappingQualityRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff

rm -f $PREFIX_LOCAL.vcf
