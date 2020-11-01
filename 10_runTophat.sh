#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=368G   # Memory per core, use --mem= for memory per node
#SBATCH --time=10-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --job-name=Tophat
##SBATCH --partition=ksu-gen-highmem.q

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
mkdir ${here}/tophatOut
cd ${here}/tophatOut

#build index file for genome
bowtie2-build $genome ${here}/tophatOut/genome_db
#align RNAseq data to genome
module load Bowtie2
$tophat2 -p 16 -o ${here}/tophatOut -p1 genome_db $rnaSeq
echo "tophat2 -o ${here}/tophatOut -p1 genome_db $rnaSeq"

FILTERBAM=
#filtering raw alignments 
samtools sort -n ${here}/tophatOut/accepted_hits.bam > ${here}/tophatOut/accepted_hits.sf.bam
FILTERBAM=/homes/mgheffel/software/augustus-3.3.3/bin/filterBam
$FILTERBAM --uniq --in ${here}/tophatOut/accepted_hits.s.bam --out ${here}/tophatOut/accepted_hits.sf.bam
samtools view -H ${here}/tophatOut/accepted_hits.sf.bam > ${here}/tophatOut/header.txt

#creating introns
samtools sort ${here}/tophatOut/accepted_hits.sf.bam both.ssf
BAM2HINTS=/homes/mgheffel/software/augustus-3.3.3/bin/bam2hints 
$BAM2HINTS --intronsonly --in=both.ssf.bam --out=${here}/tophatOut/hints.gff
