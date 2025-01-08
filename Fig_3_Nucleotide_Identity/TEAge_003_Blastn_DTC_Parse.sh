#!/bin/bash
#Program:
#SBATCH -J process
#SBATCH -n 2 # number of cores
#SBATCH --mem=10GB
#SBATCH --array=1-400%50

set -e -u -o pipefail

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)

# Relaxed
cat Blastn_Relaxed_${name}.tab |grep -v "^#" |cut -f1,2,3,4,6,7 > /Data/Fig_3_Nucleotide_Identity/small_table/Blastn_Relaxed_DTC_Small_Table/Blastn_Relaxed_${name}.small.tab

echo "# FINISH ${name} in ${HOSTNAME}"
