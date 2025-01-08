#!/bin/bash 
#Program:
#SBATCH -J harvest
#SBATCH -n 1
#SBATCH --mem=100GB
#SBATCH --array=1-6%6

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list4)
Name=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $1}')
Genome=$(echo $INPUT | awk '{FS="\t";OFS="\t"}{print $2}')

#cwd: ~/TE/03_TE_annotation/02_find_intact_LTR/02_harvest

mkdir ${Name}
cd ${Name}

#1 LTRharvest
#################################################

module load genometools/1.6.1

# build index
## LTRharvest needs the corresponding tables: tis, suf, lcp, des, ssp and sds.
## dna: processing DNA-sequences.
### ref: http://genometools.org/documents/ltrharvest.pdf

gt suffixerator \
	-db ${Genome} \
	-indexname ${Name}.gtindex \
	-tis -suf -lcp -des -ssp -sds \
	-dna
echo "Index Done"

#################################################
LTR_fa=${Name}".fa"
LTR_inner_fa=${Name}".inner.fa"
LTR_gff=${Name}".gff3"
LTR_out=${Name}".scn"

gt ltrharvest -index ${Name}.gtindex \
-seed 20 \
-minlenltr 100 \
-maxlenltr 7000 \
-similar 80 \
-mintsd 4 \
-maxtsd 6 \
-motif TGCA \
-motifmis 1 \
-vic 10 \
-out ${LTR_fa} \
-outinner ${LTR_inner_fa} \
-seqids yes \
-gff3 ${LTR_gff} > ${LTR_out}

echo "DONE"
