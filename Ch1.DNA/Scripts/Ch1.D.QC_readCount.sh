## Parameters for Apocrita job submission
    # Will simply be commented out if not running on apocrita

#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-1
#$ -tc 1


####             Ch1_M.            ####
# Written by: James B


# Run from script directory 


####        Environment Preparation         ####

## Set directories
BAMDIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.2_data/bams

## Loading sample list to create array
sample_list=/data/SBCS-eizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Metadata/JB1R_SampleList.txt

## Pulling out sample IDs
    # cat =                 # sed =
ID=$( cat "$sample_list" | sed -n ${SGE_TASK_ID}p )



#### Running Script ####

## Load samtools module
module load samtools/1.19.2

## Run samtools view to see number of reads
total_reads=$(samtools view -c ${BAMDIR}/${ID}_deduplicated.sorted_by_name.bam)
echo "Total reads for $SAMPLE: $total_reads"

## Unload module
module unload samtools/1.19.2