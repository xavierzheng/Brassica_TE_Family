#!/bin/bash

# AIM: generate the input fasta with shorter header for FIMO

set -e -u -o pipefail

echo "#1. transform fasta header to chrom::start-end format"

for name in $(cat ../name.list)
do
	cat ../IntackTE_${name}.fa | awk 'BEGIN{FS="::";OFS=""}{if($1~/^>/){print ">", $2}else{print $0}}' > IntackTE_${name}_edited.fa
done

echo "# 2. generate name table to inter-change the result"

for name in $(cat ../name.list)
do
	cat ../IntackTE_${name}.fa | awk 'BEGIN{FS="::";OFS=""}{if($1~/^>/){print $0, "\t", $2}}' |sed 's/^>//g' > IDTable_${name}.tsv
done

echo "# FINISH"
