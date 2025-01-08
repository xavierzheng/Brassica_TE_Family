#!/bin/bash 
#Program: ID_intactLTR
#SBATCH -J ID_intactLTR
#SBATCH -n 1
#SBATCH --mem=50GB
#SBATCH --array=1-12%6

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list4)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')
GenomeID=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $3}')

#cwd: /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/04_ID_intactLTR
mkdir ${Name}/
cd ${Name}/

#1 parse retriever gff and generate ID
for PasslistGFF3 in $(ls -1 /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}*.list.gff3); do
	if [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.mod.pass.list.gff3" ]; then
		
		for LTR in Copia Gypsy unknown; do
		
		less ${PasslistGFF3} | \
		grep -P "\trepeat_region" | grep -P "LTR/${LTR};" \
		> ${PasslistGFF3}.${LTR} 
		
		LINE=$( less ${PasslistGFF3}.${LTR} | wc -l )
		for ((i=1;i<=$LINE;i++)); do
		
		ltrID=$( echo $LTR | cut -c1 )
		less ${PasslistGFF3}.${LTR} | \
		sort -k1,1d -k4,4n | \
		awk -v x="${ltrID}" -v GenomeID="${GenomeID}" \
		'BEGIN{FS="\t";OFS="\t"}{print $1,$4-1,$5,GenomeID"i"x sprintf("%04d",NR)}' | \
		sort -k1,1d -k2,2n > ${Name}.${LTR}.bed

		done
		rm ${PasslistGFF3}.${LTR} 
		done
		
	elif [ -e "/nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/03_retriever/${Name}/${Fasta}.pass.list.gff3" ]; then
		
		for LTR in Copia Gypsy unknown; do
		
		less ${PasslistGFF3} | \
		grep -P "\trepeat_region" | grep -P "LTR/${LTR};" \
		> ${PasslistGFF3}.${LTR} 
		
		LINE=$( less ${PasslistGFF3}.${LTR} | wc -l )
		for ((i=1;i<=$LINE;i++)); do
		
		ltrID=$( echo $LTR | cut -c1 )
		less ${PasslistGFF3}.${LTR} | \
		sort -k1,1d -k4,4n | \
		awk -v x="${ltrID}" -v GenomeID="${GenomeID}" \
		'BEGIN{FS="\t";OFS="\t"}{print $1,$4-1,$5,GenomeID"i"x sprintf("%04d",NR)}' | \
		sort -k1,1d -k2,2n > ${Name}.${LTR}.bed

		done
		rm ${PasslistGFF3}.${LTR} 
		done
		
	else
		echo "ERROR: no pass.list.gff3 for intact Gypsy/Copia"
	fi
done

#2 generate intact LTRs bed file 
cat ${Name}.Copia.bed ${Name}.Gypsy.bed ${Name}.unknown.bed | sort -k1,1d -k2,2n > ${Name}.N.bed

#3 make gff ID file
less ${Name}.N.bed | \
awk '{FS="\t"; OFS="\t"}{print $1":"$2+1".."$3,$4 }' > ${Name}.N.gffID

echo "DONE"
