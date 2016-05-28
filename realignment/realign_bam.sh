#!/bin/bash
# Check input parameters
if [ $# -ne 3 ]
  then
	echo "Realigns a BAM file around indels"
    echo "USAGE: realign_bam.sh INPUT_BAM OUTPUT_BAM REFERENCE"
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
PREFIX_LOCAL=`basename $1`
echo $INPUT_BAM
OUTPUT_BAM=$2
echo $OUTPUT_BAM
REFERENCE=$3
echo $REFERENCE

# Preprocessing pipeline
java -Xmx1g -jar $GATK -T RealignerTargetCreator -R $REFERENCE -I $INPUT_BAM -o $PREFIX_LOCAL.intervals --filter_mismatching_base_and_quals 
java -Xmx4g -Djava.io.tmpdir=. -jar $GATK -T IndelRealigner -I $INPUT_BAM -R $REFERENCE -targetIntervals $PREFIX_LOCAL.intervals -o $OUTPUT_BAM --filter_mismatching_base_and_quals --consensusDeterminationModel USE_READS

rm -f $PREFIX_LOCAL.intervals
