#!/bin/bash

set -e -u -o pipefail

module purge
module load bedtools2/2.30.0
module load seqtk
module load seqkit


for i in $(cat name.list)
do
	cat /Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactFullInfo.gff3 |awk -v PLANT_ID=${i} 'BEGIN{FS="\t";OFS="\t"}{if($3=="helitron"){split($9, a, ";");gsub("ID=", "", a[1]); gsub("Name=", "", a[2]); print $1, $4-1, $5, PLANT_ID";"$1";"a[1]";"a[2]";DNA/Helitron", ".", $7}}' >> helitron.bed

done

for i in $(cat name.list)
do
bedtools getfasta  -s -nameOnly -fi /Data/Fig_4_TE_Feature/Genome_Fasta/${i}.genome.fa -bed helitron.bed >> helitron.fa
done

sed -i 's/([+-])//g' helitron.fa

# GC
seqtk comp helitron.fa | awk 'BEGIN{FS="\t";OFS="\t"; print "FASTA_HEAD", "LEN", "GC"}{print $1, $2, ($4+$5)/$2}' > helitron.gc

seqkit locate -d -p CG --bed helitron.fa > helitron.methyl_CG.bed
seqkit locate -d -p CHH --bed helitron.fa > helitron.methyl_CHH.bed
seqkit locate -d -p CHG --bed helitron.fa > helitron.methyl_CHG.bed

