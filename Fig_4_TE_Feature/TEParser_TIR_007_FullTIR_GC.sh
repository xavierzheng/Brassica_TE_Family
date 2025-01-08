#!/bin/bash

module purge
module load seqtk
module load seqkit

for TE in DTA DTC DTH DTM DTT
do
	seqtk comp FullLength.${TE}.all.fa | awk 'BEGIN{FS="\t";OFS="\t"; print "FASTA_HEAD", "LEN", "GC"}{print $1, $2, ($4+$5)/$2}' > FullLength_TIR_${TE}.gc

	seqkit locate -d -p CG --bed FullLength.${TE}.all.fa > TIR_${TE}.methyl_CG.bed
	seqkit locate -d -p CHH --bed FullLength.${TE}.all.fa > TIR_${TE}.methyl_CHH.bed
	seqkit locate -d -p CHG --bed FullLength.${TE}.all.fa > TIR_${TE}.methyl_CHG.bed
done
