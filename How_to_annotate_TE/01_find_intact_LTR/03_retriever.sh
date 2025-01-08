#!/bin/bash 
#Program:
#SBATCH -J retriever
#SBATCH -n 12
#SBATCH --mem=50GB
#SBATCH --array=1-6%2

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list7)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')

#cwd: /projects/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever

mkdir ${Name}
cd ${Name}

# cat finder and harvest results
cat /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/01_finder/${Name}/${Fasta}.finder.combine.scn \
/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/02_harvest/${Name}/${Name}.scn > ${Name}.rawLTR.scn

# run retriever
module load LTR_retriever/2.9.0
LTR_retriever \
	-genome ${Genome} \
	-inharvest ${Name}.rawLTR.scn \
	-threads 12 \
	-blastplus /software/shared/apps/blast/2.8.1/bin \
	-hmmer /software/shared/apps/hmmer/3.3/bin \
	-u 1.5e-8 
	
#rm *size.list *homo* *scn

cd ../
echo "DONE"


