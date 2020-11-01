#!/bin/bash
#
#       originally written by Matthew Heffel
#       11/05/2019
#
#converts fasta files to fastq files from the firsctory $HERE/raw_trimmedContigs
#fasta directroy can be changed in pyscripts/fastaTofastq.py on line 14 and 19
#
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --job-name=fastaToQ
#SBATCH --output=6C_contigsToFastq.out
#SBATCH --error=6C_contigsToFastq.err

usage() { echo "Usage: bash $0 -h data_directory(ex: /bulk/dmarth027/run_022118)" 1>&2; exit 1;}

while getopts ":h:" o; do
        case "${o}" in
                        h)
                                HERE=${OPTARG}
                                if [ ! -d $HERE ]; then echo "Directory not found $HERE"; exit 1;fi
                                HERE=$(realpath -s $HERE)
                                ;;
                        *)
                                usage
                                ;;
        esac
done
if [ -z $HERE ]; then usage; fi

module load Python
source /homes/mgheffel/virtualenvs/mypy/bin/activate
python /homes/mgheffel/algaePipeline/pyscripts/fastaToFastq.py $HERE

touch ${HERE}/6C.done
