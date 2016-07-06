#!/bin/bash

# Check input parameters
if [ $# -ne 3 ]
	then
		echo "Realigns a BAM file around indels"
    	echo "USAGE: realign_bam.sh INPUT_BAM OUTPUT_BAM REFERENCE"
		exit 1
fi

# Configuration
SCRIPT_LOCAL_2=$(readlink -f "$BASH_SOURCE")
BASEDIR_LOCAL_2=$(dirname "$SCRIPT")
source $BASEDIR_LOCAL_2/../config/config.sh

INPUT_BAM_LOCAL_2=$1
PREFIX_LOCAL_2=`basename $1 .bam`
# echo $INPUT_BAM_LOCAL_2
OUTPUT_BAM_LOCAL_2=$2
# echo $OUTPUT_BAM_LOCAL_2
REFERENCE_LOCAL_2=$3
# echo $REFERENCE_LOCAL_2
OUTPUT_DIR_LOCAL_2=$(dirname "$OUTPUT_BAM_LOCAL_2")

# Preprocessing pipeline
echo "GATK RealignerTargetCreator"
echo
java -Xmx1g -jar $GATK -T RealignerTargetCreator -R $REFERENCE_LOCAL_2 -I $INPUT_BAM_LOCAL_2 -o $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.intervals --filter_mismatching_base_and_quals
echo
echo "GATK IndelRealigner"  
echo
java -Xmx4g -Djava.io.tmpdir=. -jar $GATK -T IndelRealigner -I $INPUT_BAM_LOCAL_2 -R $REFERENCE_LOCAL_2 -targetIntervals $OUTPUT_DIR_LOCAL_2/$PREFIX_LOCAL_2.intervals -o $OUTPUT_BAM_LOCAL_2 --filter_mismatching_base_and_quals --consensusDeterminationModel USE_READS
echo
