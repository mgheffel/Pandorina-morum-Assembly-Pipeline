#!/bin/bash
#       originally written by Matthew Heffel
#	11/29/2019
#run python script to gather reference species from blast results and place them in a directory called bacterialBlastSequences
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --job-name=getBacSeqs

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

mkdir ${HERE}/bacterialBlastSequences


module load Python
#hard coded location to python virtual environment
source /homes/mgheffel/virtualenvs/mypy/bin/activate
#hard coded location to python script. Should be in $pipelineDirectory/pyscripts/
python /homes/mgheffel/algaePipeline/pyscripts/pullBacBlastSeqs.py $HERE
