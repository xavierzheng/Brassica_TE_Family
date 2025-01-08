#!/bin/bash
#Program: mcl
#SBATCH -J mcl_TEsuper9

#mkdir 03_mcl/

# 1. prepare fasta: 02_strucTEfa/ACgenome12.strucTE.fa
rm -f 02_strucTEfa/ACgenome12.strucTE.fa
for CODE in DA DB DN DP DQ DR DU IC IG IH II IJ ; do
cat 02_strucTEfa/${CODE}.strucTE.fa >> 02_strucTEfa/ACgenome12.strucTE.fa
done

# 2. prepare mcl input fasta by TE superfamily
for CLASS in $(cat 01_TEbed/*.strucTE.bed | cut -f4 | cut -f2 -d "|"| sort | uniq); do
TE=$(echo $CLASS | cut -f2 -d "/")

less 02_strucTEfa/ACgenome12.strucTE.fa | \
grep "^>" | grep ${CLASS} | \
sed 's/^>//' > 02_strucTEfa/ACgenome12.struc_${TE}.id.list.tmp

perl /nfs/project1/lab_script/group_tools/extract_sequence.pl 02_strucTEfa/ACgenome12.strucTE.fa 02_strucTEfa/ACgenome12.struc_${TE}.id.list.tmp > 03_mcl/ACgenome12.struc_${TE}.fa

rm 02_strucTEfa/ACgenome12.struc_${TE}.id.list.tmp
done

