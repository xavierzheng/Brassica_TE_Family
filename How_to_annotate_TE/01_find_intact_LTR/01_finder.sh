#!/bin/bash 
#Program:finder
#SBATCH -J finder
#SBATCH -n 10
#SBATCH --mem=50GB
#SBATCH --array=1-5%5

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list11)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')

#cwd: ~/TE/03_TE_annotation/02_find_intact_LTR/01_finder

mkdir ${Name}
cd ${Name}

#2 LTR_finder
#################################################
# run LTR_finder

## harvest_out: Output LTRharvest format if specified. Default: LTR_FINDER table format.
## size: size to split the genome sequence. Default 5,000,000 (bp).
## time: Default: 1500 (sec). Suggestion: 300 for -size 1000000.
### ref: https://github.com/oushujun/LTR_FINDER_parallel

module load LTR_FINDER_parallel/1.1

LTR_FINDER_parallel -seq ${Genome} -threads 8 -harvest_out -size 72500000 -time 21750

echo "DONE"
