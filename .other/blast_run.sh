#!/bin/bash
#Refer to  blastn.sh for running with default parameters
#STEP 1: 	Copy this file to your directory
#STEP 2:	Change [dollar]1 to the path of your fasta file
#STEP 3:	Change parameters in the last command as you see fit
#STEP 4:	Run:
#			qsub 'blastn_file_default.sh file.fasta'

#SBATCH --mem=150G   # Memory per core, use --mem= for memory per node
#SBATCH --time=10-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=blastn

usage() { echo "sbatch $0 to_blast.fasta blast_output.blast" 1>&2; exit 1; }

#if (($# < 2));then usage;fi

QUERY=$1
OUT=${2-"${1}.blastn"}
BLASTN=/homes/bioinfo/bioinfo_software/ncbi-blast-2.6.0+/bin/blastn
export BLASTDB="/bulk/bioinfo/bioinfo_software/blastdb/main"

if [ ! -e $QUERY ]; then echo "File does not exist $QUERY"; usage; fi

$BLASTN -db nt_all -num_threads $SLURM_NTASKS  -penalty -2 -reward 1 -word_size 28 -outfmt "6 qseqid evalue sacc sseqid length pident bitscore stitle" -query "$QUERY" -out "$OUT" -max_target_seqs 1

echo "$BLASTN -db nt_all -num_threads $SLURM_NTASKS  -penalty -2 -reward 1 -word_size 28 -outfmt "6 qseqid evalue sacc sseqid length pident bitscore stitle" -query "$QUERY" -out "$OUT" -max_target_seqs 1"
