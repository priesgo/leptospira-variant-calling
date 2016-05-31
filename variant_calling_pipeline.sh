#!/bin/bash
# Check input parameters
if [ $# -ne 2 ]
  then
	echo "Run the whole variant calling pipeline and outputs the union and intersection VCFs"
    echo "USAGE: variant_calling_pipeline.sh <INPUT_BAM> <OUTPUT_FOLDER> <REFERENCE>"
    exit 1
fi

#mkdir logs
#> logs/$PREFIX.log
#2> logs/$PREFIX.log

# Parameters
SCRIPT=$(readlink -f "$BASH_SOURCE")
VARIANT_CALLING_BASEDIR=$(dirname "$SCRIPT")
INPUT_BAM=$1
PREFIX=`basename $1 .bam`
echo $INPUT_BAM
OUTPUT_FOLDER=$2
echo "Output folder: $OUTPUT_FOLDER"
REFERENCE=$3
echo $REFERENCE
REFERENCE_BASENAME=`basename $REFERENCE`

declare -A SNPEFF_REFERENCES=( ["Lb.Hardjo.L550.fasta"]="Leptospira_borgpetersenii_serovar_Hardjo_bovis_L550_uid58507" ["Lb.Hardjo.JB197.fasta"]="Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509")
SNPEFF_REFERENCE="${SNPEFF_REFERENCES[$REFERENCE_BASENAME]}"

# BAM preprocessing
source $VARIANT_CALLING_BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_FOLDER/$PREFIX.preprocessed.bam

# Realignment around indels
source $VARIANT_CALLING_BASEDIR/realignment/realign_bam.sh $OUTPUT_FOLDER/$PREFIX.preprocessed.bam $OUTPUT_FOLDER/$PREFIX.realigned.bam $REFERENCE

# Calls variants with Haplotype Caller
source $VARIANT_CALLING_BASEDIR/variant_calling/haplotype_caller.sh $OUTPUT_FOLDER/$PREFIX.realigned.bam $OUTPUT_FOLDER/$PREFIX.hc.raw.vcf $REFERENCE

# Filters false positive variants
source $VARIANT_CALLING_BASEDIR/variant_filtering/variant_filtering.sh $OUTPUT_FOLDER/$PREFIX.hc.raw.vcf $OUTPUT_FOLDER/$PREFIX.hc.filtered.vcf $REFERENCE

# Calls variants with samtools
source $VARIANT_CALLING_BASEDIR/variant_calling/samtools_pileup.sh $OUTPUT_FOLDER/$PREFIX.realigned.bam $OUTPUT_FOLDER/$PREFIX.st.filtered.vcf $REFERENCE

# Filters false positive variants
# Variants are already filtered in the samtools script
#source $BASEDIR/variant_filtering/variant_filtering.sh $PREFIX.st.raw.vcf $PREFIX.st.filtered.vcf $REFERENCE

# Calls variants with Unified Genotyper
source $VARIANT_CALLING_BASEDIR/variant_calling/unified_genotyper.sh $OUTPUT_FOLDER/$PREFIX.realigned.bam $OUTPUT_FOLDER/$PREFIX.ug.raw.vcf $REFERENCE

# Filters false positive variants
source $VARIANT_CALLING_BASEDIR/variant_filtering/variant_filtering.sh $OUTPUT_FOLDER/$PREFIX.ug.raw.vcf $OUTPUT_FOLDER/$PREFIX.ug.filtered.vcf $REFERENCE

# combines variants from all callers
source $VARIANT_CALLING_BASEDIR/combine_variants/combine_variants.sh $OUTPUT_FOLDER/$PREFIX.hc.filtered.vcf $OUTPUT_FOLDER/$PREFIX.ug.filtered.vcf $OUTPUT_FOLDER/$PREFIX.st.filtered.vcf $OUTPUT_FOLDER/$PREFIX.union.vcf $OUTPUT_FOLDER/$PREFIX.intersection.vcf $REFERENCE

# annotation
source $VARIANT_CALLING_BASEDIR/annotation/annotation.sh $OUTPUT_FOLDER/$PREFIX.intersection.vcf $OUTPUT_FOLDER/$PREFIX.intersection.annotated.vcf $SNPEFF_REFERENCE $REFERENCE

source $VARIANT_CALLING_BASEDIR/annotation/annotation.sh $OUTPUT_FOLDER/$PREFIX.union.vcf $OUTPUT_FOLDER/$PREFIX.union.annotated.vcf $SNPEFF_REFERENCE $REFERENCE
