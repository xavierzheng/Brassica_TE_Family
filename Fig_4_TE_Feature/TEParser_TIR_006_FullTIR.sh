#!/bin/bash
# program:
#       bedtools
#SBATCH -J bedtools
#SBATCH -n 4
#SBATCH --mem=4G


module purge
module load bedtools2/2.30.0

for i in $(cat name.list)
do
	for TE in DTA DTC DTH DTM DTT
	do
		cat /Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactSingleLine.sort.gff3 | grep "Classification=TIR/${TE}" | awk -v PLANT_ID=${i} 'BEGIN{FS="\t";OFS="\t"}{split($9, a, ";"); gsub("ID=", "", a[1]); gsub("Name=", "", a[2]); gsub("Classification=", "", a[3]); print $1, $4-1, $5, PLANT_ID";"$1";"a[1]";"a[2]";"a[3], ".", $7}' > FullLength.${TE}.${i}.bed
	done
done

echo "# GET Fasta========"

for i in $(cat name.list)
do
	for TE in DTA DTC DTH DTM DTT
	do
		bedtools getfasta -s -nameOnly -fi /Data/Fig_4_TE_Feature/Genome_Fasta/${i}.genome.fa -bed FullLength.${TE}.${i}.bed >> FullLength.${TE}.${i}.fa
	done
done

echo "# Put together========"

for TE in DTA DTC DTH DTM DTT
do
	cat FullLength.${TE}.D*.bed FullLength.${TE}.I*.bed >> FullLength.${TE}.all.bed
	rm FullLength.${TE}.D*.bed FullLength.${TE}.I*.bed
	cat FullLength.${TE}.D*.fa FullLength.${TE}.I*.fa >> FullLength.${TE}.all.fa
	rm FullLength.${TE}.D*.fa FullLength.${TE}.I*.fa
done

