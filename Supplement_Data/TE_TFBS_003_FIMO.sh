#!/bin/bash
#Program:
#SBATCH -J fimo
#SBATCH -n 1 # number of cores
#SBATCH --mem=10GB
#SBATCH --array=1-50%8

set -e -u -o pipefail

INPUT=$1

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list.split_${INPUT})

WD="/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN/TF_FIMO"

motif="/nfs/project_ssd/project3/pxzhe/JASPAR/20231006080751_JASPAR2022_combined_matrices_4417_meme.txt"

# split query fasta in order speed up 
# for i in $(cat name.list); do mkdir SPLIT_${i} ;cd SPLIT_${i}; ln -s ../../IntackTE_${i}.fa . ; perl /homes/pxzhe/bioscript/equalSplitFasta.pl IntackTE_${i}.fa 50 . ; cd ../ ;done
# for i in $(cat name.list); do cd SPLIT_${i} ; ls IntackTE_${i}.fa_* > ../name.list.split_${i}; cd ../ ; done

module purge
module load meme/5.5.2

fimo --oc ${WD}/FIMO_OUT_${name} \
    --verbosity 2 \
    --bfile bg_markov_${INPUT}.model \
    --text \
    ${motif} SPLIT_${INPUT}/${name} > FIMO_OUT_${name}

echo "# FINISH ${name} in ${HOSTNAME}"
