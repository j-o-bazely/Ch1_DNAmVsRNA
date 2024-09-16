##### Ch1.D.4_GTFtoBED.sh   ####

#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y

## This script takes in a GTF file and converts it to a BED12 file format, used in genome anntation

##IMPORTANT: Run script from the 'Scripts' directory in which it is stored

##          Environment Preparation         ##
# Move into cleanPHD to set WORKDIR path from which all other required directories can be accessed downstream
cd ../../..

## Set paths
WORKDIR=$(pwd)
DATADIR="{$WORKDIR}Ch1_DataStorage/"

echo "Working directory: $WORKDIR"
echo "Data directory: $DATADIR"

##              Script to Run               ##

## Move into directory containing Gff3 file
cd {DATADIR}/Ch1_inputData/Ch1_ChangGenome

echo "Now in $(pwd)" 

# Run gff3 to BED conversion
$TRANSDECODER_PATH/gff3_file_to_bed.pl \
Chang23LoggerheadAnnotation.gff \
> Chang23LoggerheadAnnotation.bed12
