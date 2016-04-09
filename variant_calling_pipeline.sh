# Check input parameters
if [ $# -ne 2 ]
  then
    echo "USAGE: variant_calling_pipeline.sh $INPUT_BAM $REFERENCE"
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
source $BASEDIR/haplotype_caller/haplotype_caller.sh $PREFIX.realigned.bam $PREFIX.raw.vcf $REFERENCE

# Filters false positive variants
source $BASEDIR/variant_filtering/variant_filtering.sh $PREFIX.raw.vcf $PREFIX.filtered.vcf $REFERENCE

# annotation
source $BASEDIR/annotation/annotation.sh $PREFIX.filtered.vcf $PREFIX.annotated.vcf $SNPEFF_REFERENCE

