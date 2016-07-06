#!/bin/bash
# Check input parameters
if [ $# -ne 1 ]
  then
	echo "Creates all required indices (i.e.: fai and dict). Input file must have .fasta extension"
	echo "USAGE: prepare_reference.sh INPUT_FASTA"
	exit 1
fi

# Configuration
SCRIPT_REF=$(readlink -f "$BASH_SOURCE")
BASEDIR_REF=$(dirname "$SCRIPT_REF")
source $BASEDIR_REF/../config/config.sh

INPUT_REF=$1
# echo $INPUT_REF
PREFIX_REF=`basename $1 .fasta`
OUTPUT_DIR_REF=$(dirname "$INPUT_REF")

echo "Preparing the reference ..."
echo

# FAI index
echo "$SAMTOOLS_HOME/samtools faidx $INPUT_REF"
$SAMTOOLS_HOME/samtools faidx $INPUT_REF

# Picard's dict index
echo "java -jar $PICARD CreateSequenceDictionary R=$INPUT_REF O=${OUTPUT_DIR_REF}/${PREFIX_REF}.dict"
java -jar $PICARD CreateSequenceDictionary R=$INPUT_REF O=${OUTPUT_DIR_REF}/${PREFIX_REF}.dict

echo
