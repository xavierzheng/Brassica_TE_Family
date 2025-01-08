#!/bin/bash
#Program:
#SBATCH -J parse_age
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH --array=1-31%16

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list1)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')
GenomeID=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $3}')

#cwd: /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/06_age

mkdir ${Name}
cd ${Name}

# 1. get age from retriever passlist
for Passlist in $(ls -1 /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}*.pass.list); do
	if [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.mod.pass.list" ]; then
		less /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.mod.pass.list | \
		cut -f1,5,6,10,12 -d "$(printf '\t')" | \
		awk '{FS="\t"}NR>1{OFS="\t"; split($1,a,":"); split($2,b,"."); split($3,c,"."); \
		if ($2 == ".." && $3 == "..") print a[1]":"a[2], $5+1, $4; \
		else print a[1]":"b[1]".."c[3], $5+1, $4}' > ${Name}.age.name

	elif [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.pass.list" ]; then
		less /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.pass.list | \
		cut -f1,5,6,10,12 -d "$(printf '\t')" | \
		awk '{FS="\t"}NR>1{OFS="\t"; split($1,a,":"); split($2,b,"."); split($3,c,"."); \
		if ($2 == ".." && $3 == "..") print a[1]":"a[2], $5+1, $4; \
		else print a[1]":"b[1]".."c[3], $5+1, $4}' > ${Name}.age.name
	else
		echo "ERROR: no pass.list for intactLTR's age"
	fi
done

# 2. change intact LTR IDs.
less ${Name}.age.name | \
awk 'BEGIN{FS="\t";OFS="\t"}{if (NR==FNR)d[$1]=$2; \
else{$1=d[$1];print $0}}' /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/04_ID_intactLTR/${Name}/${Name}.N.gffID - > tmp
mv tmp ${Name}.age.name



# 2. change name
#module load Python/3.9.4
#scp ~/script/name.pkl ./
#scp ~/script/change_LTR_name.py ./
#python3 change_LTR_name.py ${Name}.age
#rm name.pkl
#rm change_LTR_name.py

echo "DONE"
