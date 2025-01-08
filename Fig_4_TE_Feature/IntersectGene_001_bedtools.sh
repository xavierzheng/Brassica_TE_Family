#!/bin/bash

set -e -u -o pipefail

module purge
module load bedtools2/2.30.0

for i in $(cat name.list)
do

zcat /Data/Fig_4_TE_Feature/GFF_Gene_TE/GeneAnnotate_${i}.InChr.sort.gff.gz |awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="gene"){split($9, a, ";");print $1, $2, $3, $4, $5, $6, $7, $8, a[1]}}' > GeneOnly.${i}.gff

GENE_GFF="GeneOnly.${i}.gff"
TE_GFF="/Data/Fig_4_TE_Feature/GFF_Gene_TE/${i}.TE.InChr.IntactSingleLine.sort.gff3"

bedtools intersect -a ${TE_GFF} -b ${GENE_GFF} -wo > Gene_IntactTE_relation.${i}.tab

done
