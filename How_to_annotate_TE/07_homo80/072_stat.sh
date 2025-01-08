#!/bin/bash
#Program: samtools
#SBATCH -J count
#SBATCH -n 1
#SBATCH --mem=50GB
#SBATCH --array=1-5%5
module purge
module load blast/2.10.1
module load bedtools/2.29.2

Genome=$(sed -n "$SLURM_ARRAY_TASK_ID"p genome.list)
GenomeFasta=$(less /homes/kochiaying/name32.list | grep -w $Genome | cut -f1)
Name=$(less /homes/kochiaying/name32.list | grep -w $Genome | cut -f7 -d "/")
Gff=$(echo "/nfs/project1/Repeat/10_integrate/01_CACTA/"$Name"/"$Name".EDTA_scanCACTA.gff3")

cd $Name

# 1. dictionary replace
#for REPEAT in LTRCO LTRGY LTRRT ; do

#for REPEAT in LTRCO LTRGY LTRRT DTA DTH DTT DTM Helitron ; do
#Bed=$(echo "/nfs/project1/Repeat/10_integrate/07_homo80/BrassicaAC_2/"$Name"/"$Name".EDTA_scanCACTA."$REPEAT".bed")

# make a list of pass TE
#less $Name.$REPEAT.TE_family.list | \
#cut -f2 | tr ';' '\n' | awk 'NF>0' | \
#sort | uniq > $Name.$REPEAT.homo80.list.tmp

# remove homo-TEs (NOstrucTE and didnt pass 80-80-80)
#rm -f $Name.EDTA_scanCACTA.$REPEAT.homo80.bed.tmp
#for PASS in $(less $Name.$REPEAT.homo80.list.tmp); do
#less $Bed | \
#awk -v PASS=$PASS 'BEGIN{FS="\t";OFS="\t"}{split($4,a,";");\
#if(a[1]==PASS)print}' >> $Name.EDTA_scanCACTA.$REPEAT.homo80.bed.tmp
#done

# dictionary replace bed file
#less $Name.EDTA_scanCACTA.$REPEAT.homo80.bed.tmp | \
#sed 's/_INT;/;/;s/_LTR;/;/' | \
#awk 'BEGIN{FS="\t";OFS="\t"}{split($4,a,";");\
#if(NR==FNR)d[$1]=$2;\
#else{gsub(a[2],d[a[2]],$4);print}}' $Name.$REPEAT.TE_Name.list - 

# remove tmp
#rm $Name.$REPEAT.homo80.list.tmp
#rm $Name.EDTA_scanCACTA.$REPEAT.homo80.bed.tmp

# print gff
#done | sort -k1,1 -k2,2n -k3,3n | \
#awk 'BEGIN{FS="\t";OFS="\t"}\
#{split($4,a,";");\
#if(a[1] ~/^TE_homo/){METH="homology"}else{METH="structural"};\
#print $1,"homo80",$4,$2+1,$3,$5,$6,".","ID="a[1]";Name="a[2]";Classification="a[3]";Methodology="METH}' | \
#sort -k1,1 -k4,4n -k5,5n > $Name.homo80.gff.tmp


# replace TE family name of CACTA
#less /nfs/project1/Repeat/10_integrate/06_clusterTE/CC/04-parse_id/TE_family.list | \
#grep "^CACTA" | awk 'BEGIN{FS="\t";OFS="\t"}\
#{split($1,TEfam,"_");split($2,a,"(");split(a[1],NAME,":");split(NAME[4],NUM,"-");\
#print NAME[3]":"NUM[1]+1"-"NUM[2],TEfam[1]}' > CACTA.dict.tmp

#less $Gff | \
#awk BEGIN'{FS="\t";OFS="\t"}\
#{if ($3 == "intact_CACTA") print}' | \
#awk 'BEGIN{FS="\t";OFS="\t"}{if(NR==FNR){d[$1]=$2}\
#else{split($9,a,"Name=");split(a[2],NAME,";");\
#gsub(NAME[1],d[NAME[1]],$9);print}}' CACTA.dict.tmp - > $Name.DTC.gff.replaced.tmp

#less $Name.DTC.gff.replaced.tmp | \
#awk 'BEGIN{FS="\t";OFS="\t"}\
#{split($9,a,"ID=");split(a[2],ID,";");\
#split($9,b,"Classification=");split(b[2],CLASS,";");\
#split($9,c,"Name=");split(c[2],NAME,";");\
#print $1,$2,ID[1]";"NAME[1]";"CLASS[1],$4,$5,$6,$7,$8,$9}' > $Name.DTC.gff.tmp


# cat gff with CACTA
#cat $Name.homo80.gff.tmp $Name.DTC.gff.tmp | \
#grep -v "Name=;Classification=" | sort -k1,1 -k4,4n -k5,5n > $Name.EDTA_scanCACTA.homo80.gff

#rm $Name.homo80.gff.tmp $Name.DTC.gff.tmp 
#rm CACTA.dict.tmp $Name.DTC.gff.replaced.tmp

# 2. stat table
Out="$Name.EDTA_scanCACTA.homo80.gff" 
sh /nfs/project1/Repeat/10_integrate/00_stat_TEclass_gff.sh ${Out} ${GenomeFasta}.fai > ../stat/stat_${Name}.table


echo "DONE"

