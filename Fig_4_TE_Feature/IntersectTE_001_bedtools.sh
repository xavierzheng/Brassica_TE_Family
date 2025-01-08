#!/bin/bash

set -e -u -o pipefail

module purge
module load bedtools2/2.30.0

for i in $(cat name.list)
do

	INTACT_TE="/Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactSingleLine.sort.gff3"
	TOTAL_TE="/Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactSingleLineAndHomoTE.sort.gff3"

	# remove the intersect results which contain themselves 

	bedtools intersect -a ${INTACT_TE} -b ${TOTAL_TE} -wo | awk 'BEGIN{FS="\t";OFS="\t"}{if($4==$13 && $5==$14 && $9==$18){}else{print $0}}' - > IntactTE_TotalTE_relation.${i}.tab
done

