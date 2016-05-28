# Check input parameters
if [ $# -ne 3 ]
  then
    echo "Preprocess BAM and realigns it around indels."
    echo "USAGE: prepare_bam.sh INPUT_BAM OUTPUT_FOLDER REFERENCE"
    exit 1
fi

#mkdir logs
#> logs/$PREFIX.log
#2> logs/$PREFIX.log

# Parameters
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
INPUT_BAM=$1
PREFIX=`basename $1`
echo "Input BAM: $INPUT_BAM"
OUTPUT_FOLDER=$2
echo "Output folder: $OUTPUT_FOLDER"
REFERENCE=$3
echo "Reference: $REFERENCE"
REFERENCE_BASENAME=`basename $REFERENCE`


# BAM preprocessing
source $BASEDIR/preprocessing/preprocess_bam.sh $INPUT_BAM $OUTPUT_FOLDER/$PREFIX.preprocessed.bam

# Realignment around indels
source $BASEDIR/realignment/realign_bam.sh $OUTPUT_FOLDER/$PREFIX.preprocessed.bam $OUTPUT_FOLDER/$PREFIX.realigned.bam $REFERENCE
