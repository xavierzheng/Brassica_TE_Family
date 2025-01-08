#!/bin/bash
# Program:
#       run blastx
#SBATCH -J blastx
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH --array=1-1%1

set -e -u -o pipefail

INPUT=$1

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list.split_${INPUT})

module load blast/2.13.0

QUERY="/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN/TF_FIMO/SPLIT_${INPUT}/${name}"
BlastDATABASE="/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/Compare_REXdb/Viridiplantae_v3.0_ALL_protein-domains_repet_formated.fa"
OutPutFile="BLASTX_REXdb_${name}.out"

blastx \
    -query $QUERY \
    -db $BlastDATABASE \
    -out $OutPutFile \
    -query_gencode 1 \
    -outfmt "7 qseqid sseqid qlen slen sstrand pident length mismatch gapopen gaps qframe sframe qstart qend sstart send evalue bitscore" \
    -max_target_seqs 500 \
    -num_threads 1 \
    -evalue 1e-5 \
    -soft_masking false

echo "${name} in $HOSTNAME"
