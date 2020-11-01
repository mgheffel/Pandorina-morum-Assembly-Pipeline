#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=150G   # Memory per core, use --mem= for memory per node
#SBATCH --time=10-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --job-name=Augustus
module load bzip2
module load ncurses

export PATH=$PATH:/homes/bioinfo/bioinfo_software/bowtie2-2.3.4.3/
module load SAMtools 
tophat2=/homes/mgheffel/software/tophat-2.1.1.Linux_x86_64/tophat2

usage() { echo "Usage: bash $0 -g genome.fa -r RNAseq.fa -h $PWD" 1>&2; exit 1; }

while getopts ":g:r:h:" o; do
	case "${o}" in
		g)		
			genome=${OPTARG};;
		r)
			rnaSeq=${OPTARG};;
		h)
			here=${OPTARG};;
		*)
			usage;;
	esac
done
if [ -z $genome ] || [ -z $rnaSeq ]; then usage; fi
#mkdir ${here}/tophatOut
cd ${here}/tophatOut
module load AUGUSTUS
augustus --species=chlamydomonas --codingseq=on --alternatives-from-evidence=true --hintsfile=${here}/tophatOut/hints.gff --extrinsicCfgFile=/homes/mgheffel/software/augustus-3.3.3/config/extrinsic/extrinsic.cfg --allow_hinted_splicesites=atac ${here}/tophatOut/genome_db.fa > augustus.out
