#!/bin/bash
#Program:
#SBATCH -J EDTA
#SBATCH -n 21 
#SBATCH --mem=80GB
#SBATCH --array=1-30%10

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list1)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')

#Name="Test3"
#Genome="/nfs/project1/Brassica/Genome_reID/Boleracea_CC/2021_1892/1892_reID.fa"
#Fasta="1892_reID.fa"

set -euo pipefail

module purge
module load singularity/3.5.2

#################################################################################

mkdir ${Name}
cd ${Name}

scp ${Genome} ./

singularity exec -B /nfs/project1,/usr/lib/locale/:/usr/lib/locale/ \
	/software/shared/apps/tools/docker_singularity/EDTA/1.9.5/EDTA.sif EDTA.pl --genome ${Fasta} \
	--overwrite 0 \
	--sensitive 1 \
	--anno 1 \
	--step all \
	--evaluate 1 \
	--threads 7

echo "DONE"
