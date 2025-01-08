#!/bin/bash
#Program: getfasta
#SBATCH -J getfasta
#SBATCH -n 1
#SBATCH --mem=50GB
#SBATCH -x chopin
#SBATCH --array=1-12%12

module purge
module load bedtools/2.29.2
module load htslib/1.13

#mkdir 01_TEbed
#mkdir 02_strucTEfa

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)
CODE=$(echo $INPUT | awk 'BEGIN{FS=" "}{print $2}')
Genome=$(echo $INPUT | awk 'BEGIN{FS=" "}{print $1}')
Name=$(echo $INPUT | cut -f7 -d "/")
GFF_strucTE=$(\
echo "/nfs/project1/Repeat/10_integrate/03_helitron/"${Name}"/"${Name}".strucTE.gff3")

# 0. decompress genome
#bgzip -d ${Genome}.gz

# 1. strucTE bed
less ${GFF_strucTE}| \
awk 'BEGIN{FS="\t";OFS="\t"}\
{if ($3=="Copia_LTR_retrotransposon" || $3=="Gypsy_LTR_retrotransposon" || $3=="LTR_retrotransposon" || $3=="Mutator_TIR_transposon"|| $3=="PIF_Harbinger_TIR_transposon" || $3=="Tc1_Mariner_TIR_transposon" || $3=="hAT_TIR_transposon" || $3=="helitron" || $3=="intact_CACTA") print}' | \
awk 'BEGIN{FS="\t";OFS="\t"}\
{split($9,a,"Classification=");split(a[2],b,";");\
print $1,$4-1,$5,$1":"$4"-"$5"|"b[1],$6,$7}' | \
sort -k1,1 -k2,2n > 01_TEbed/${CODE}.strucTE.bed


# 2. get fasta of strucTE
bedtools getfasta -nameOnly \
        -fi ${Genome} \
        -bed 01_TEbed/${CODE}.strucTE.bed \
        -fo 02_strucTEfa/${CODE}.strucTE.fa

# 3.make a good fasta format
less 02_strucTEfa/${CODE}.strucTE.fa | \
grep "^>" | sed 's/^>//' > 02_strucTEfa/${CODE}.strucTE.id.tmp  

perl /nfs/project1/lab_script/group_tools/extract_sequence.pl 02_strucTEfa/${CODE}.strucTE.fa 02_strucTEfa/${CODE}.strucTE.id.tmp > 02_strucTEfa/${CODE}.strucTE.tmp
mv 02_strucTEfa/${CODE}.strucTE.tmp 02_strucTEfa/${CODE}.strucTE.fa
rm 02_strucTEfa/${CODE}.strucTE.id.tmp

echo ${CODE}"-DONE"
