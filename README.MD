This bioinformatic pipeline was used to assemble the v1 nuclear genome of Pandorina morum. The link to the manuscript will be added at publication.
This pipeline was built and run on Kansas State University's High-Performance Computing cluster, Beocat. Beocat is a
Beowult cluster of CentOS Linux servers coordinated by the Slurm job scheduling system.

#---NOTE---# many scripts contain hard coded paths which will need to be changed for external use.
Three subdirectories exist; pyscripts/ contains python scripts used in main numbered submission scripts and 
.other/ contains a script for blast submission used by 2_blast.sh and a side script for calculating the n50 of a fasta file.
abyss/ contains scripts for running ABySS as steps 1 and 7 of the pipeline

Steps 1 and 7 are not labeled by their script title and involve running the paired-end sequence assembler ABySS,
the process of which is described here. Abyss can be run without using the provided scripts here.
a) chunk_pair5.pl
-Seperates the r1 and r2 files into n files of roughly equal size (we used 16)
Usage chunk_pair5.pl --read1 FILE1.fastq --read2 FILE2.fastq --nfiles n --library_name LIBRARY_NAME --library_type (Single|Paired|Mate-paired) --logfile LOGFILE.txt

b) abyss-master3.pl
Creates and runs ABySS scripts in parallel for different k-mer sizes
Usage: abyss-master3.pl --data_file logfile_fromChunkPair
-more options at bottom of script

Steps 2-6, and 8-12 can be completed by simply running the scripts in order with their desired parameters.

2_blast.sh 
-Runs nucleotide blast on files in input directory with blast options:
"db nt_all -num_threads $SLURM_NTASKS  -penalty -2 -reward 1 -word_size 28 -outfmt "6 qseqid evalue sacc sseqid length pident bitscore stitle""
Usage: 2_blast.sh -i directory_of_blastable_seqeuences [-a blast_configuration*] [-c job_configuration]
Files will output into `$input | rev | cut -f 2- -d '_' | rev)_blast`

3_runPullBacSeqs.sh
-Runs python script to gather reference species from blast results and place them in a directory called bacterialBlastSequences
Usage: 3_runPullBacSeqs.sh -h $PWD
Script expects blast output files to be in a directory $PWD/scaffolds_blast/ and outputs reference sequences 
into a directory called $PWD/bacterialBlastSequences
The species to be obtained may be changed from bacteria by altering the regex in the script pullBacBlastSeqs.py

4_remapLargeMemPerf.sh
-Runs Bowtie2 remaping a forward and reverse read to a reference sequence with perfect matches only
Usage: 4_remapLargeMemPerf.sh -f reference.fasta -f forwardRead.fastq -r reverseReads.fastq -o outputDir
Output of Bowtie2 job will be placed in $outputDir/

5_runRemoveReads.sh
-Runs python script that removes reads from forward and reverse read file if the read appears in the input sam (alignment) file
Usage: 5_runRemoveReads.sh -h $PWD -s remappingFile.sam -f forwardRead.fastq -r reverseRead.fastq
Output goes to r1name_noBac.fastq and r2name_noBac.fastq

6_runFixDupeReads.sh (optional)
-Runs a python script that adds extensions to duplicate read names to keep both reads but prevent suplicate ID
This script exists because at this point in my dataset I merged different sequencing runs into one file.
Usage: 6_runFixDupeReads.sh -i input.fastq
Output is input-noDupe.fastq

7 (ABySS stage)

8_contigsToFastq
-expects $PWD/raw_trimmedContigs/ to be a directory with fasta files (directory name can be changed in python script)
-Converts all fasta files in the directory to fastq files with perfect quality score
Usage: 8_contigsToFastq -h $PWD
Output files go to $PWD/raw_fastqContigs

9A_canuCorrect.sh (optional)
-Runs correct stage of Canu
Usage 9A_canuCorrect.sh -r singleMoleculeLongReads.fastq -p prefix -d assemblyDir -s genomeSize
Output goes to $assemblyDir

9B_canuTrim.sh (optional)
-Runs trim stage of Canu
Usage 9B_canuTrim.sh -r singleMoleculeLongReads.fastq -p prefix -d assemblyDir -s genomeSize
Output goes to $assemblyDir

9C_canuAssemble.sh
-Runs assemble stage of Canu
Usage 9C_canuAssemble.sh.sh -r singleMoleculeLongReads.fastq -p prefix -d assemblyDir -s genomeSize
Output goes to $assemblyDir

10_runTophat.sh
-Runs tophat to map RNAseq reads to a reference (in this case the newly assembled genome)
Usage: 10_runTophat.sh -g genome.fasta -r rnaSeq.fastq -h $PWD
Output goes to $PWD/tophatOut

11_makeAugustusHints.sh
Converts tophat output into a hints.gff file for use with Augustus
Usage: 11_makeAugustusHints.sh -h $PWD
expects input to be in $PWD/tophatOut/
Output is $PWD/tophatOut/himts.gff

12_augustus.sh
Runs Augustus using the genome, rnaSeq data and the hints.gff file
Usage: 12_augustus.sh -g genome.fasta -r rnaSeq.fastq -h $PWD
expects hints file to in $PWD/tophatOut/hints.gff
Output is augustus.out

Downstream analysis was done using orthoMCL, PosiGene, and a Linear Regression Model to view homologous group expansion/contraction
No scripts werre used for OrthoMCL, it was run using all default parameters and generated homologous gene groups.
PosiGene was used only to generate latest common ancestor(LCA) distances using one outgroup.