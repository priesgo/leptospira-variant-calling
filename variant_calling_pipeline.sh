# Check input parameters
if [ $# -ne 2 ]
  then
	echo "Run the whole variant calling pipeline and outputs the union and intersection VCFs"
    echo "USAGE: variant_calling_pipeline.sh INPUT_BAM REFERENCE"
    exit 1
fi

#mkdir logs
#> logs/$PREFIX.log
#2> logs/$PREFIX.log

# Parameters
BASEDIR=$(dirname "$0")
INPUT_BAM=$1
PREFIX=`basename $1`
echo $INPUT_BAM
REFERENCE=$2
echo $REFERENCE
REFERENCE_BASENAME=`basename $REFERENCE`

declare -A SNPEFF_REFERENCES=( ["Lb.Hardjo.L550.fasta"]="Leptospira_borgpetersenii_serovar_Hardjo_bovis_L550_uid58507" ["Lb.Hardjo.JB197.fasta"]="Leptospira_borgpetersenii_serovar_Hardjo_bovis_JB197_uid58509")
SNPEFF_REFERENCE="${SNPEFF_REFERENCES[$REFERENCE_BASENAME]}"

# BAM preprocessing
source $BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $PREFIX.preprocessed.bam

# Realignment around indels
source $BASEDIR/realignment/realign_bam.sh $PREFIX.preprocessed.bam $PREFIX.realigned.bam $REFERENCE

# Calls variants with Haplotype Caller
source $BASEDIR/haplotype_caller/haplotype_caller.sh $PREFIX.realigned.bam $PREFIX.hc.raw.vcf $REFERENCE

# Filters false positive variants
source $BASEDIR/variant_filtering/variant_filtering.sh $PREFIX.hc.raw.vcf $PREFIX.hc.filtered.vcf $REFERENCE

# Calls variants with samtools
source $BASEDIR/samtools_pileup/samtools_pileup.sh $PREFIX.realigned.bam $PREFIX.st.filtered.vcf $REFERENCE

# Filters false positive variants
# Variants are already filtered in the samtools script
#source $BASEDIR/variant_filtering/variant_filtering.sh $PREFIX.st.raw.vcf $PREFIX.st.filtered.vcf $REFERENCE

# Calls variants with Unified Genotyper
source $BASEDIR/unified_genotyper/unified_genotyper.sh $PREFIX.realigned.bam $PREFIX.ug.raw.vcf $REFERENCE

# Filters false positive variants
source $BASEDIR/variant_filtering/variant_filtering.sh $PREFIX.ug.raw.vcf $PREFIX.ug.filtered.vcf $REFERENCE

# combines variants from all callers
source $BASEDIR/combine_variants/combine_variants.sh $PREFIX.hc.filtered.vcf  $PREFIX.ug.filtered.vcf $PREFIX.st.filtered.vcf $PREFIX.union.vcf $PREFIX.intersection.vcf $REFERENCE

# annotation
source $BASEDIR/annotation/annotation.sh $PREFIX.intersection.vcf $PREFIX.intersection.annotated.vcf $SNPEFF_REFERENCE $REFERENCE

source $BASEDIR/annotation/annotation.sh $PREFIX.union.vcf $PREFIX.union.annotated.vcf $SNPEFF_REFERENCE $REFERENCE
