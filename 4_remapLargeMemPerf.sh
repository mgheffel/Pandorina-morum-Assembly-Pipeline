#!/bin/bash -l


#SBATCH --mem=256G
#SBATCH --time=2-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
##SBATCH --partition=killable.q,batch.q
#SBATCH --job-name=remap_bowtie

module load SAMtools
module load bzip2
module load ncurses

PATH=$PATH:/homes/bioinfo/bioinfo_software/bowtie2-2.3.3.1:/homes/bioinfo/bioinfo_software/samtools-1.7/bin

usage() { echo "sbatch $0 -d reference -f forward_read -r reverse_read -o output_directory" &1>2;exit 1;}

while getopts ":d:f:r:o:" o;do
	case "${o}" in
		d)
			ref=${OPTARG}
			if [ ! -e $ref ]; then echo "File does not exist $ref";usage;fi
			ref=$(realpath $ref)
			;;
		f)	
			r1=${OPTARG}
			if [ ! -e $r1 ]; then echo "File does not exist $r1";usage; fi
			r1=$(realpath $r1)
			;;
		r)	
			r2=${OPTARG}
			if [ ! -e $r2 ]; then echo "File does not exist $r2";usage; fi
			r2=$(realpath $r2)
			;;
		o)	
			outdir=${OPTARG}
			;;
		*)
			echo "Parameter not supported : ${OPTARG}"
			usage
			;;

	esac
done
if [ -z $ref ] || [ -z $r1 ] || [ -z $r2 ] || [ -z $outdir ];then echo "Requeired parameter missing";usage;fi

#create flag to schow script is in progress
touch ${outdir}${ref}.hold

dir=$(dirname $r1)
mkdir $outdir
outdir=$(realpath $outdir)
sample=$outdir/$(basename $r1 | cut -f 1 -d '.')

mkdir -p ${sample}_index
bowtie2-build -f $ref ${sample}_index/$(basename $sample)

bowtie2 --score-min 'C,0,-1' -x ${sample}_index/$(basename $sample) -1 $r1 -2 $r2 --un-conc ${sample}_unmapped.fastq -p 4 > $sample.sam

samtools view -bS $sample.sam > $sample.bam
samtools sort $sample.bam -o $sample.sorted.bam
samtools index $sample.sorted.bam

#remove flag
rm -rf ${outdir}${ref}.hold
