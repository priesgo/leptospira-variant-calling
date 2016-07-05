#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
	echo "Runs the whole variant calling pipeline and outputs the union and intersection VCFs"
	echo "USAGE: variant_calling_pipeline.sh <INPUT_BAM> <OUTPUT_DIR> <REFERENCE>"
	exit 1
fi

# Input parameters
SCRIPT=$(readlink -f "$BASH_SOURCE")
VARIANT_CALLING_BASEDIR=$(dirname "$SCRIPT")

#echo $VARIANT_CALLING_BASEDIR
INPUT_BAM=$1
echo "Input BAM file: $INPUT_BAM"
PREFIX=`basename $1 .bam`
OUTPUT_DIR=$2
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir $OUTPUT_DIR
fi
echo "Output directory: $OUTPUT_DIR"
REFERENCE=$3
echo "Reference genome: $REFERENCE"
REFERENCE_BASENAME=`basename $REFERENCE`
echo

# Prepare the reference
source $VARIANT_CALLING_BASEDIR/reference/prepare_reference.sh $REFERENCE

# Preparing input BAM
source $VARIANT_CALLING_BASEDIR/prepare_bam.sh $INPUT_BAM $OUTPUT_DIR $REFERENCE

# Variant calling with samtools
source $VARIANT_CALLING_BASEDIR/variant_calling/samtools_pileup.sh $OUTPUT_DIR/$INPUT_BAM.realigned.bam $OUTPUT_DIR/samtools/$PREFIX.raw.vcf $REFERENCE

# Variant calling with GATK
source $VARIANT_CALLING_BASEDIR/variant_calling/unified_genotyper.sh $OUTPUT_DIR/$INPUT_BAM.realigned.bam $OUTPUT_DIR/unified_genotyper/$PREFIX.raw.vcf $REFERENCE
source $VARIANT_CALLING_BASEDIR/variant_calling/haplotype_caller.sh $OUTPUT_DIR/$INPUT_BAM.realigned.bam $OUTPUT_DIR/haplotype_caller/$PREFIX.raw.vcf $REFERENCE

# Filter variants using GATK
source $VARIANT_CALLING_BASEDIR/variant_filtering/variant_filtering_samtools.sh $OUTPUT_DIR/samtools/$PREFIX.raw.vcf $OUTPUT_DIR/samtools/$PREFIX.filtered.vcf $REFERENCE
source $VARIANT_CALLING_BASEDIR/variant_filtering/variant_filtering_gatk.sh $OUTPUT_DIR/haplotype_caller/$PREFIX.raw.vcf $OUTPUT_DIR/haplotype_caller/$PREFIX.filtered.vcf $REFERENCE
source $VARIANT_CALLING_BASEDIR/variant_filtering/variant_filtering_gatk.sh $OUTPUT_DIR/unified_genotyper/$PREFIX.raw.vcf $OUTPUT_DIR/unified_genotyper/$PREFIX.filtered.vcf $REFERENCE

# Combine variants from the three variant callers
source $VARIANT_CALLING_BASEDIR/combine_variants/combine_variants.sh $OUTPUT_DIR/haplotype_caller/$PREFIX.filtered.vcf $OUTPUT_DIR/unified_genotyper/$PREFIX.filtered.vcf $OUTPUT_DIR/samtools/$PREFIX.filtered.vcf variants/$PREFIX.union.vcf variants/$PREFIX.intersection.vcf $REFERENCE
