## Parameters for Apocrita job submission
    # Will simply be commented out if not running on apocrita

#!/bin/bash
#$ -pe smp 3
#$ -l h_vmem=10G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-17
#$ -tc 1


####             Ch1_M.            ####
# Written by: James B


# Run from script directory 


####        Environment Preparation         ####

## Set directories
BAMDIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.ReMap_data/Ch1.ReMap.2_data/bams

## Loading sample list to create array
sample_list=/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Metadata/Ch1.SampleList.txt

## Pulling out sample IDs
    # cat =                 # sed =
ID=$( cat "$sample_list" | sed -n ${SGE_TASK_ID}p )



#### Running Script ####

## Load samtools module
module load samtools/1.19.2

## Run samtools view to see number of reads
total_reads=$(samtools view -c ${BAMDIR}/${ID}/${ID}_deduplicated.sorted_by_name.bam)
echo "Total reads for $ID: $total_reads"


## Unload module
module unload samtools/1.19.2
