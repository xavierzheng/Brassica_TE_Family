#!/bin/bash 
#Program: scan_to_gff
#SBATCH -J scan_to_gff
#SBATCH -n 1
#SBATCH --array=1-5%5
#SBATCH --mem=50GB

INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list5)
Name=$(echo $INPUT | awk '{print $1}')
Genome=$(echo $INPUT | awk '{print $2}')

module purge
module load scan_for_matches/20211222

mkdir ${Name}/
cd ${Name}/

# 1. scan CACTA
scan_for_matches /homes/kochiaying/CACTA/pattern2 \
<${Genome}> ${Name}_scanCACTA.fa

# 2. get gff
less ${Name}_scanCACTA.fa | \
grep '^>' | sed 's/^>//;s/:/,/;s/\[//;s/\]//' | \
awk 'BEGIN{FS=",";OFS="\t"}{print $1,"scan_for_matches","intact_CACTA",$2,$3,"0","+",".",
"Name="$1":"$2"-"$3";Classification=TIR/DTC;Method=structural"}' | \
awk 'BEGIN{FS="\t";OFS="\t"}\
{print $1,$2,$3,$4,$5,$6,$7,$8,"ID=CACTA_"sprintf("%04d",NR)";"$9}\
{print $1,$2,"CACTA_lTSD",$4,$4+2,$6,$7,$8,"ID=CACTA_lTSD"sprintf("%04d",NR)";"$9}\
{print $1,$2,"CACTA_lTIR",$4+3,$4+14,$6,$7,$8,"ID=CACTA_lTIR"sprintf("%04d",NR)";"$9}\
{print $1,$2,"CACTA_internal",$4+15,$5-15,$6,$7,$8,"ID=CACTA_internal"sprintf("%04d",NR)";"$9}\
{print $1,$2,"CACTA_rTIR",$5-14,$5-3,$6,$7,$8,"ID=CACTA_rTIR"sprintf("%04d",NR)";"$9}\
{print $1,$2,"CACTA_rTSD",$5-2,$5,$6,$7,$8,"ID=CACTA_rTSD_"sprintf("%04d",NR)";"$9}' | \
sort -k1,1 -k4,4n > ${Name}.intactCACTA.gff3

echo "DONE"





