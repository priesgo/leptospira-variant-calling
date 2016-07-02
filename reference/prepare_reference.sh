#!/bin/bash
# Check input parameters
if [ $# -ne 1 ]
  then
	echo "Creates all required indices (i.e.: fai and dict). Input file must have .fasta extension"
    echo "USAGE: prepare_reference.sh INPUT_FASTA"
    exit 1
fi

#Configuration
SCRIPT=$(readlink -f "$BASH_SOURCE")
BASEDIR=$(dirname "$SCRIPT")
source $BASEDIR/../config/config.sh

# Copies input BAM into local folder
INPUT_REFERENCE=$1
echo $INPUT_REFERENCE
PREFIX_LOCAL=`basename $1 .fasta`
OUTPUT_DIR=$(dirname "$INPUT_REFERENCE")

# FAI index
$SAMTOOLS_HOME/samtools faidx $INPUT_REFERENCE

# Picard's dict index
java -jar $PICARD CreateSequenceDictionary R=$INPUT_REFERENCE O=${OUTPUT_DIR}/${PREFIX_LOCAL}.dict

