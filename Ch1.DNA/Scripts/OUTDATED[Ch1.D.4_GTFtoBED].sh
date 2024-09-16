##### Ch1.D.4_GTFtoBED.sh   ####

#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y

## This script takes in a GTF file and converts it to a BED12 file format, used in genome anntation

## Load module
module load bedops/2.4.35

## Set paths
genomeDir=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1_inputData/Ch1_ChangGenome/

## Move into genome folder
cd $genomeDir

## Run GTF to BED12 conversion
gtf2bed < genomic.gtf > genomic.bed12

## Unload the bedops module
module unload bedops/2.4.35