#!/bin/bash
#Program:
#SBATCH -J markov
#SBATCH -n 4 # number of cores
#SBATCH --mem=10GB
#SBATCH --array=2-9%4

set -e -u -o pipefail

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../name.list)

WD="/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN/TF_FIMO"

module purge
module load meme/5.5.2

cd ${WD}

fasta-get-markov -dna -m 0 IntackTE_${name}_edited.fa bg_markov_${name}.model

echo "# FINISH ${name} in ${HOSTNAME}"
