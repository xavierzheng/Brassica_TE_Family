#!/bin/bash 
#Program:
#SBATCH -J sum_intact
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH --array=1-12%6

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list4)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')
Fasta=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}'| awk 'BEGIN{FS="/"}{print $8}')
GenomeID=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $3}')

#cwd: ~/TE/03_TE_annotation/02_find_intact_LTR/03_retriever
cd ${Name}/

# 1. copy steps 01 and 02 result here
scp /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/01_finder/${Name}/${Fasta}.finder.combine.gff3 ./${Name}.finder.gff3
scp /nfs/project1/Repeat/03_TE_annotation/01_find_intact_LTR/02_harvest/${Name}/${Name}.gff3 ./${Name}.harvest.gff3

# 2. calculate 
F1=$( \
less ${Name}.finder.gff3 | grep -P "LTR_FINDER_parallel\trepeat_region" | wc -l \
)

H2=$( \
less ${Name}.harvest.gff3 | grep -P "LTRharvest\trepeat_region" | wc -l \
)

R3_c=$( \
less ${Name}.N.bed | cut -f 4 | awk 'BEGIN{FS="i"}{if ($2 ~/^C/) print $0}' | wc -l \
)

R3_g=$( \
less ${Name}.N.bed | cut -f 4 | awk 'BEGIN{FS="i"}{if ($2 ~/^G/) print $0}' | wc -l \
)

R3_u=$( \
less ${Name}.N.bed | cut -f 4 | awk 'BEGIN{FS="i"}{if ($2 ~/^u/) print $0}' | wc -l \
)

R3=$( \
less ${Name}.N.bed | wc -l \
)

echo -e ${Name}"\t"${Fasta}"\t"${F1}"\t"${H2}"\t"${R3_c}"\t"${R3_g}"\t"${R3_u}"\t"${R3} >> ../stat_intact.txt

# 3. remove redundant
#rm ${Name}.finder.gff3
#rm ${Name}.harvest.gff3

echo "DONE"



