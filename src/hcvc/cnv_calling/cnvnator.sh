#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Runs CNVnator"
    echo "USAGE: cnvnator.sh INPUT_BAM OUTPUT_FOLDER REFERENCE"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Copies input BAM into local folder
#cp $1 .
#cp $1.bai .
INPUT_BAM=$1
echo $INPUT_BAM
PREFIX_LOCAL=`basename $1 .bam`
OUTPUT_FOLDER=$2
echo $OUTPUT_FOLDER
REFERENCE=$3
echo $REFERENCE


# Prepare reference genome: split reference by chromosome and copy in destination folder
$GENOMETOOLS_HOME/gt splitfasta -splitdesc $OUTPUT_FOLDER $REFERENCE

# Runs CNVnator pipeline
$CNVNATOR_HOME/cnvnator -root $PREFIX_LOCAL.NC_008508.root -genome $REFERENCE -chrom NC_008508 -tree $INPUT_BAM 
$CNVNATOR_HOME/cnvnator -root $PREFIX_LOCAL.NC_008508.root -genome $REFERENCE -chrom NC_008508 -his 300
$CNVNATOR_HOME/cnvnator -root $PREFIX_LOCAL.NC_008508.root -genome $REFERENCE -chrom NC_008508 -stat 300
$CNVNATOR_HOME/cnvnator -root $PREFIX_LOCAL.NC_008508.root -genome $REFERENCE -chrom NC_008508 -partition 300
$CNVNATOR_HOME/cnvnator -root $PREFIX_LOCAL.NC_008508.root -genome $REFERENCE -chrom NC_008508 -call 300 > $PREFIX_LOCAL.NC_008508.cnvnator


# Runs BAQ with PrintReads as HaplotypeCaller does not support it on the fly as UnifiedGenotyper does
echo "GATK PrintReads to calculate BAQ"
java -jar $GATK -T PrintReads -R $REFERENCE -I $INPUT_BAM -baq RECALCULATE -o $OUTPUT_DIR/$PREFIX_LOCAL.baq.bam

# Haplotype caller variant calling pipeline
echo "GATK HaplotypeCaller"
java -jar $GATK -T HaplotypeCaller -R $REFERENCE -I $OUTPUT_DIR/$PREFIX_LOCAL.baq.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 --min_base_quality_score 13 --downsampling_type NONE -ploidy 1 -nda -allowNonUniqueKmersInRef -bamout $INPUT_BAM.hc_reassembly.bam -o $OUTPUT_VCF --annotateNDA --annotation BaseQualityRankSumTest --annotation ClippingRankSumTest --annotation Coverage --annotation FisherStrand --annotation GCContent --annotation HomopolymerRun --annotation LikelihoodRankSumTest --annotation NBaseCount --annotation QualByDepth --annotation RMSMappingQuality --annotation StrandOddsRatio --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample --annotation DepthPerSampleHC --annotation StrandAlleleCountsBySample --annotation StrandBiasBySample --excludeAnnotation HaplotypeScore --excludeAnnotation InbreedingCoeff
# Use this parameters to create all haplotypes -forceActive -disableOptimizations

rm -f $OUTPUT_DIR/$PREFIX_LOCAL.baq.ba*
