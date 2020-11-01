#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=100G   # Memory per core, use --mem= for memory per node
#SBATCH --time=10-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --job-name=Augustus
CWD="/homes/mgheffel/algaePipeline"
module load bzip2
module load ncurses

export PATH=$PATH:/homes/bioinfo/bioinfo_software/bowtie2-2.3.4.3/
module load SAMtools 
tophat2=/homes/mgheffel/software/tophat-2.1.1.Linux_x86_64/tophat2

usage() { echo "Usage: bash $0 -h $PWD" 1>&2; exit 1; }

while getopts ":h:" o; do
	case "${o}" in
		h)
			here=${OPTARG};;
		*)
			usage;;
	esac
done
if [ -z $here ]; then usage; fi
#mkdir ${here}/tophatOut
cd ${here}/tophatOut
module load Bowtie2
echo "loaded bowtie"
echo "sorting accepted_hits"
samtools sort -n ${here}/tophatOut/accepted_hits.bam > ${here}/tophatOut/accepted_hits.s.bam
echo "filtering sorted accepted_hits"
FILTERBAM=/homes/mgheffel/software/augustus-3.3.3/bin/filterBam
$FILTERBAM --uniq --in ${here}/tophatOut/accepted_hits.s.bam --out ${here}/tophatOut/accepted_hits.sf.bam
echo "samtools view"
samtools view -H ${here}/tophatOut/accepted_hits.sf.bam > ${here}/tophatOut/header.txt

#creating introns
echo "creating introns"
samtools sort ${here}/tophatOut/accepted_hits.sf.bam > ${here}/tophatOut/both.ssf.bam
BAM2HINTS=/homes/mgheffel/software/augustus-3.3.3/bin/bam2hints 
$BAM2HINTS --intronsonly --in=${here}/tophatOut/both.ssf.bam --out=${here}/tophatOut/hints.gff
