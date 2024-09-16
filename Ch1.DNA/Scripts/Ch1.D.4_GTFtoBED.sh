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
WORKDIR="{$pwd}/"
DATADIR="{$WORKDIR}Ch1_DataStorage/"

echo "Working directory: {$WORKDIR}"
echo "Data directory: {$DATADIR}"

## Load module

