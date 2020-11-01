#!/bin/bash
#       originally written by Matthew Heffel
#       11/04/2019
#ensures reads do not have duplicate IDs after merging read sets

#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --job-name=rmvReads

usage() { echo "Usage: bash $0 -i input_fastqFile" 1>&2; exit 1;}

while getopts ":i:" o; do
        case "${o}" in
			i)
				inFile=${OPTARG}
				;;
                        *)
                                usage
                                ;;
        esac
done
if [ -z $inFile ] ; then usage; fi

module load Python
#hard coded path to python virtual environment
source /homes/mgheffel/virtualenvs/mypy/bin/activate
#hard coded path to python script
python /homes/mgheffel/algaePipeline/pyscripts/fixDupeReads.py $inFile
