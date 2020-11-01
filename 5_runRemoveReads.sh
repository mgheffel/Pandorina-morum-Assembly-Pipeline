#!/bin/bash
#       originally written by Matthew Heffel
#       10/24/2019

#runs python script to remove reads from a forward and reverse file that are matched in a .bam file to a species declared in the python script (Bacteria superkingdom)

#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --job-name=rmvReads

usage() { echo "Usage: bash $0 -h $PWD -s remappingFile.sam -f forwardRead.fastq -r reverseRead.fastq" 1>&2; exit 1;}

while getopts ":h:s:f:r:" o; do
        case "${o}" in
                        h)
                                HERE=${OPTARG}
                                if [ ! -d $HERE ]; then echo "Directory not found $HERE"; exit 1;fi
                                HERE=$(realpath -s $HERE)
                                ;;
			s)
				samFile=${OPTARG}
				;;
			f)
				r1File=${OPTARG}
				;;
			r)
				r2File=${OPTARG}
				;;
                        *)
                                usage
                                ;;
        esac
done
if [ -z $HERE ] || [ -z $samFile ]; then usage; fi

module load Python
#hard coded pah to python virtual environment
source /homes/mgheffel/virtualenvs/mypy/bin/activate
#hard coded path to script
python /homes/mgheffel/algaePipeline/pyscripts/removeMappedReads.py $samFile $r1File $r2File
