#!/bin/bash

## obtain long_terminal_repeat from GFF
## remove Parent=repeat_region_28;
for i in $(cat name.list); do cat /Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactFullInfo.gff3 |awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="long_terminal_repeat"){print $0}}' |sed 's/Parent=repeat_region_[0-9]\{1,\};//g' |awk 'BEGIN{FS="\t";OFS="\t"}{split($9, a, ";"); print $1, $2, $3, $4, $5, $6, $7, $8, a[1]";"a[2]";"a[3]}' > ${i}.LongTerminalRepeat.gff3 ; done

## bedtools getfasta
module purge
module load bedtools2/2.30.0

for i in $(cat name.list); do cat ${i}.LongTerminalRepeat.gff3 |awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $4-1, $5, $9, ".", $7}' > ${i}.LongTerminalRepeat.bed ;done

for i in $(cat name.list); do bedtools getfasta -fi /Data/Fig_4_TE_Feature/Genome_Fasta/${i}.genome.fa -bed ${i}.LongTerminalRepeat.bed -name -s > ${i}.LongTerminalRepeat.fa ;done

## estimate GC
module load seqtk/1.3
for i in $(cat name.list); do seqtk comp ${i}.LongTerminalRepeat.fa |awk 'BEGIN{FS="\t";OFS="\t"}{print $1, ($4+$5)/$2}' > ${i}.LongTerminalRepeat.gc ; done

## estimate CG CHH CHG
module load seqkit/v1.3
for i in $(cat name.list); do seqkit locate -d -p CG --bed ${i}.LongTerminalRepeat.fa > ${i}.LongTerminalRepeat.methyl_CG.bed ; done 
for i in $(cat name.list); do seqkit locate -d -p CHH --bed ${i}.LongTerminalRepeat.fa > ${i}.LongTerminalRepeat.methyl_CHH.bed ; done
for i in $(cat name.list); do seqkit locate -d -p CHG --bed ${i}.LongTerminalRepeat.fa > ${i}.LongTerminalRepeat.methyl_CHG.bed ; done

## SEQ LEN
for i in $(cat name.list); do seqtk comp ${i}.LongTerminalRepeat.fa |awk 'BEGIN{FS="\t";OFS="\t"}{print $1, $2}' > ${i}.LongTerminalRepeat.length ; done

