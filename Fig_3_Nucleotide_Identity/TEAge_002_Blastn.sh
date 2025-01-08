#!/bin/bash
#Program:
#SBATCH -J blastn
#SBATCH -n 8 # number of cores
#SBATCH --mem=50GB
#SBATCH --array=1-8%8

set -e -u -o pipefail

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)

module load blast/2.13.0

QUERY="IntackTE_${name}.fa"
BlastDATABASE="IntackTE_${name}.fa"


OutPutFile="Blastn_Relaxed_${name}.tab" # v3


# ver3. use more relaxed parameters
blastn \
	-query $QUERY \
	-db $BlastDATABASE \
	-out $OutPutFile \
	-task blastn \
	-word_size 7 \
	-reward 1 \
	-penalty -1 \
	-gapopen 0 \
	-gapextend 2 \
	-outfmt "7 qseqid sseqid qlen slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads 8 \
	-evalue 1e-5 \
	-max_target_seqs 5000000 \
	-soft_masking false \
	-dust no

# -word_size 7 is blastn-short
# -reward 1
# -penalty -1
# a ratio of about one (1/-1) is best for sequences that are 75% conserved [1].
# ref: States DJ, Gish W, and Altschul SF (1991) METHODS: A companion to Methods in Enzymology 3:66-70.
# ref: https://blast.ncbi.nlm.nih.gov/doc/blast-topics/blastsearchparams.html#reward-and-penalty-for-nucleotide-programs

echo "# FINISH ${name} ${HOSTNAME}"

