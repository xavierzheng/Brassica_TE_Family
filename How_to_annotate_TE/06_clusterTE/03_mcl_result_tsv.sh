#!/bin/bash

# 4. result tsv
cd 03_mcl/

rm -f mclblastline.tsv
echo -e "TE\tInflation point\tTotal group\tn>=10 group\tn>=50 group\tn>=80 group\tSingleton\tMember" > mclblastline.tsv

for TE in $(cat ../TE.list); do
for i in 15 20 25 30 35 40 45 50 55 60 ; do
File="I_${i}/dump.out.ACgenome12.struc_${TE}.blast_out.I${i}"

Total=$(less ${File} | wc -l )
n10=$(less ${File} | awk 'BEGIN{FS="\t"}{if (NF>=10) print NF}' | wc -l )
n50=$(less ${File} | awk 'BEGIN{FS="\t"}{if (NF>=50) print NF}' | wc -l )
n80=$(less ${File} | awk 'BEGIN{FS="\t"}{if (NF>=80) print NF}' | wc -l )
Sing=$(less ${File} | awk 'BEGIN{FS="\t"}{if(NF==1)print}' | wc -l )
Memb=$(less ${File} | tr "\n" "\t"  |  sed 's/[\t]\+$//' |awk 'BEGIN{FS="\t";OFS="\t"}{print NF}' )
echo -e ${TE}"\t"${i}"\t"${Total}"\t"${n10}"\t"${n50}"\t"${n80}"\t"${Sing}"\t"${Memb} >> mclblastline.tsv

done
done
