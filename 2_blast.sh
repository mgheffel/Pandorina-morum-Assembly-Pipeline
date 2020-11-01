#!/bin/bash:
#$ -cwd'
#
#	
#
#       Updated v1.2 by Matthew G. Heffel mgheffel@ksu.edu
#       9/25/2020
#
#       Updated v1.1 by Reza Mazloom rmazloom@ksu.edu
#       10/23/2017
#
#       originally written by Majed Alsadhan
#
#

CWD="/homes/mgheffel/algaePipeline"
DIR="$CWD/.other"
script=$DIR"/blast_run.sh"
#PERL=$DIR"/get_first_hits.pl"
#FINISH=$DIR"/perl_batch_run.sh"
LOGS="$CWD/logs/"
chain="$DIR/chain.sh"
JOBS=""

usage() { echo "Usage: bash $0 -i directory_of_blastable_seqeuences [-a blast_configuration*] [-c job_configuration]" 1>&2; exit 1; }

while getopts ":i:a:c:b:" o; do
	case "${o}" in
		i)		
			input=${OPTARG}
			if [ ! -d $input ]; then echo "Invalid directory: $input";fi
			;;
		a)
			blast_args=${OPTARG}
			;;
		c)
			config=${OPTARG}
			;;
		b)
			debugDir=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done
if [ -z $input ]; then usage; fi
input=$(realpath -s $input)

mkdir $debugDir

file_list=$(ls $input/*.fasta | grep -v ".blast.sh")

RESULTS="$(echo $input | rev | cut -f 2- -d '_' | rev)_blast"
#RESULTS="${input}_blast"

rm -rf $RESULTS
mkdir $RESULTS

#cd $RESULTS

for file in $file_list
do
	file_name=$(basename $file | rev | cut -f 2- -d '.' | rev)
	echo "Blasting and submitting file: " $(basename $file)
	X=$(basename $file)"_blast"
	#JID=$(sbatch --job-name=D$X $config $script $file $RESULTS/$file_name.blast $blastargs)
	JID=$(sbatch --job-name=blastn $config $script $file $RESULTS/$file_name.blast $blastargs)
	JID=$(echo $JID | rev | cut -f 1 -d ' '|rev)
	JOBS="${JOBS}$JID,"
done
JOBS=$(echo $JOBS | sed 's/,\([^,]*\)$/ \1/')
sbatch $config --dependency $JOBS $chain $input/../7_blast.done &>/dev/null
#qsub $config -hold_jid $JOBS -N "BLAST_R" -l mem=1G,h_rt=1:00:00 -e $LOGS"BLAST_R.log" -o $LOGS"BLAST_R.log"  -pe single 1 $FINISH $PERL $RESULTS &>/dev/null
echo 
echo 
echo "Done submitting BLAST jobs" #, once the jobs are done running, another job will run to sort things out with another submission named \"BLAST_R\""
echo "You will find your blast results in $RESULTS"
echo
echo "Please wait for a couple of minutes" # until the job \"BLAST_R\" starts and finishes, then you will find the csv files"
echo 
echo 
echo "Thank you!"

