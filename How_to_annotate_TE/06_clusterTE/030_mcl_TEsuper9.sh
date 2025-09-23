#!/bin/bash
#Program: mcl
#SBATCH -J mcl_TEsuper9
#SBATCH -n 20
#SBATCH --mem=150GB
#SBATCH -x chopin
#SBATCH --array=6-6%1

module purge
module load blast/2.10.1
#cat 01_TEbed/*.strucTE.bed | cut -f4 | cut -f2 -d "|"| sort | uniq | cut -f2 -d "/" > TE.list
TE=$(sed -n "$SLURM_ARRAY_TASK_ID"p TE.list)

cd 03_mcl/

# 1. make db
makeblastdb -dbtype nucl -in ACgenome12.struc_${TE}.fa  -out ACgenome12.struc_${TE}

# 2. blastn
blastn -query ACgenome12.struc_${TE}.fa \
       -db ACgenome12.struc_${TE} \
       -task blastn \
       -out ACgenome12.struc_${TE}.align_out \
       -num_threads 20 \
       -dust no \
       -max_hsps 25000 \
       -soft_masking false

blastn -query ACgenome12.struc_${TE}.fa \
       -db ACgenome12.struc_${TE} \
       -task blastn \
       -outfmt 6 \
       -out ACgenome12.struc_${TE}.blast_out \
       -num_threads 20 \
       -dust no \
       -max_hsps 25000 \
       -soft_masking false

# 3. mclblastline
module load MCL/14-137
# --blast-m9: expect BLAST column format
# --blast-score=b: <b|e|r> (bit scores|e-values|norm bit score)
# --blast-sort=a: <a|o> (alphabetic|occurrence sorting)
# --blast-bcut=5: ignore bit scores not exceeding 5.
# --mcl-I=2.5: inflation point
# --ecut=<val> (E-value cutoff).

for j in 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 ; do
i=$(echo $j*10 | bc | sed 's/\.0$//')

mclblastline --blast-m9 --blast-score=b --blast-sort=a --blast-bcut=5 --mcl-I=${j} ACgenome12.struc_${TE}.blast_out

#mkdir I_${i}
cd I_${i}
mv ../out.ACgenome12.struc_${TE}.blast_out.I${i} ./
cd ../
done



echo ${TE}"-DONE"
