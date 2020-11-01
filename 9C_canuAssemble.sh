#!/bin/sh
#SBATCH -J CanuCorrect
#SBATCH --mem-per-cpu=256G
#SBATCH --time=3-00:00:00

usage() { echo "Usage: -r reads.fastq -p assemblyPrefix -d assemblyDir -s genomeSize" 1>&2; exit 1;}

while getopts ":r:p:d:s:" o; do
        case "${o}" in
                        r)
                                reads=${OPTARG}
                                ;;
                        p)
                                prefix=${OPTARG}
                                ;;
                        d)
                                assemblyDir=${OPTARG}
                                ;;
                        s)
                                genomeSize=${OPTARG}
                                ;;
                        *)
                                usage
                                ;;
        esac
done
if [ -z $reads ] || [ -z $prefix ] || [ -z $assemblyDir ] || [ -z $genomeSize ]; then usage; fi

echo $reads
echo $prefix
echo $assemblyDir
echo $genomeSize
hostname

#need these?
module load Java
module load gnuplot
JAVA_OPTS="-Xms12g -Xmx124g -XX:MaxPermSize=32g"
#???

/homes/bioinfo/bioinfo_software/canu-1.8/Linux-amd64/bin/canu stopOnLowCoverage=2 -p $prefix -d $assemblyDir genomeSize=$genomeSize -pacbio-raw $reads  gridOptions="--time=4:00:00"

