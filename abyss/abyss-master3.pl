#!/usr/bin/perl

# Written by Brad Olson and Dave Turner - Kansas State University

#Version 3 rewrite handles multiple libraries and store them in a data structure

use warnings;
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;

#Get the path the script is executing in, note that chunk_pair3.pl must be in this path!
my $exec_path = dirname(__FILE__); # Executable path

#Variables without default values
my($read1, $read2, $single_reads, $fasta_single_reads, $chunks);
my($prechunk_read1,$prechunk_read2,$prechunk_single,$prechunk_fasta_single); #This is a file that has chunked files in it from previous runs
my($chunk_log_1, $chunk_log_2, $chunk_log_single); #Log files for chunking
my $email; #stores email address
my $priority='-P KSU-GEN-BIOINFO';
#my $priority='-P KSU-GEN-HIGHMEM';
my %data_file;

#Variables with default values that can be reset on the command line
my $out_dir = 'abyss-assembly';
my $name = 'abyss';
my $kmer = 21;
my $kmer_increment = 5;
my $kmer_stop = 31; # Not sure if there is an upper limit for kmer anymore

#Default values for variables
my $nodes = 1;
my $cores = 16;
my $maxtime = "1";
my $total_mem = 256;


my $data_file; #File that holds all data for assembly
GetOptions ("data_file=s"		=>	\$data_file,
			"outdir=s"			=>	\$out_dir,
			"name=s"			=>	\$name,
			"kstart=i"			=>	\$kmer,
			"kstop=i"			=>	\$kmer_stop,
			"increment=i"		=>	\$kmer_increment,
			"nodes=i"			=>	\$nodes,
			"cores=i"			=>	\$cores,
			"maxtime=s"			=>	\$maxtime,
			"mem=i"				=>	\$total_mem,
			"chunks=i"			=> 	\$chunks,
			"email=s"			=>	\$email,
			"priority=s"		=>	\$priority
			) or die ("Error in command line arguments\n");

unless(defined($email)){
	die "You must provide a email address using the --email command line option";
}



#Get working directory with files
#note that pwd does not add a '/'
#add '/' manually to end
my $cwd = `pwd`;
chomp($cwd);
$cwd.='/';

#Variables set by code
#SGE has memory specified on a per core basis
my $mem_per_core = $total_mem / $cores;

#Choose the smallest version for the given kmer size
my $abyss_exec;
if( $kmer <= 32 ) {
  $abyss_exec="/homes/bjsco/local/bin/abyss-32/bin/abyss-pe";
} elsif ( $kmer <= 64 ) {
  $abyss_exec="/homes/bjsco/local/bin/abyss-64/bin/abyss-pe";
} elsif ( $kmer <= 96 ) {
  $abyss_exec="/homes/bjsco/local/bin/abyss-96/bin/abyss-pe";
} elsif ( $kmer <= 128 ) {
  $abyss_exec="/homes/bjsco/local/bin/abyss-128/bin/abyss-pe";
} elsif ( $kmer <= 256 ) {
	$abyss_exec="/homes/bjsco/local/bin/abyss-256/bin/abyss-pe";
} else {
  say "\n\nkmer size cannot be greater than 128 without recompiling\n\n";
  exit(0);
}


#Load data sets that are prepared and chunked
my $data_hash=load_data($data_file);

#Now unwind data structure to load up abyss with files to run
#Library names for each of the three types of reads Singles/Paired/Mate-pair
my @single_libraries;
my @paired_libraries;
my @mp_libraries;

for my $library_name (sort {$data_hash->{$a} <=> $data_hash->{$b}} keys %{ $data_hash }){
	my $library_type=$data_hash->{$library_name}{'Type'};
	if ($library_type eq 'Singles'){
		push (@single_libraries, $library_name);
	}
	elsif($library_type eq 'Paired'){
		push (@paired_libraries, $library_name);
	}
	elsif($library_type eq 'Mate-pair'){
		push (@mp_libraries, $library_name);
	}
	else{
		say "Warning: Library $library_name has unknown type $library_type";
	}

}

#Variables hold library names or file names
my $single_names;
my $single_files;
my $paired_names;
my $paired_files;
my $mp_names;
my $mp_files;

#Now build the abyss file sections sections
if (@single_libraries > 0){
	#$single_names="se=\'";
	foreach my $lib (@single_libraries){
		#$single_names.="$lib ";
		$single_files.="se=\'$data_hash->{$lib}{'Singles'}\' ";
	}
	#add last single quote to names
	$single_names=~s/\s$/\'/g;
}
else{ #if no single libraries
	$single_names='';
	$single_files='';
}

if (@paired_libraries > 0){
	$paired_names="lib=\'";
	foreach my $lib (@paired_libraries){
		$paired_names.="$lib ";
		$paired_files.="$lib=\'$data_hash->{$lib}{'Read1'} ";
		$paired_files.="$data_hash->{$lib}{'Read2'}\' ";
	}
	#add last single quote to names
	$paired_names=~s/\s$/\'/g;
}
else{ #if no single libraries
	$paired_names='';
	$paired_files='';
}

if (@mp_libraries > 0){
	$mp_names="mp=\'";
	foreach my $lib (@mp_libraries){
		$mp_names.="$lib ";
		$mp_files.="$lib=\'$data_hash->{$lib}{'Read1'} ";
		$mp_files.="$data_hash->{$lib}{'Read2'}\' ";
	}
	#add last single quote to names
	$mp_names=~s/\s$/\'/g;
}
else{ #if no single libraries
	$mp_names='';
	$mp_files='';
}



#Debug
# say $single_names;
# say $single_files;
# say $paired_names;
# say $paired_files;
# say $mp_names;
# say $mp_files;
####****#####
#Here we debug the output line for executing abyss
#say "abyss se=\'$single_names\' lib=\'$paired_names\'  mp=\'$mp_names\' $paired_files";
#say "Abyss line will be....";
#say "$abyss_exec $single_names $paired_names $mp_names $single_files $paired_files $mp_files"
####****####


#Now write scripts and execute them

my $kmer_dir;
my $scriptfile;

#This will make the output directory if needed
$out_dir="$cwd/$out_dir";
if( -r $out_dir ) {
   say "$out_dir already exists";
} else {
   say "Making $out_dir";
   `mkdir $out_dir`;
}

#This will create an output subdirectory, bash script, and qsub for each kmer
while ($kmer <= $kmer_stop)
{
         #Each kmer will output into its own directory

   $kmer_dir="$out_dir/kmer${kmer}";

   if( -r $kmer_dir ) {
      say "$kmer_dir already exists";
   } else {
      say "Making $kmer_dir";
      `mkdir $kmer_dir`;
   }

      #Make the individual script for this kmer

   $scriptfile = 'kmer'.$kmer.'.sh';
   say "Making script file $scriptfile";

   open(SCRIPT, "> $scriptfile");

#DDT - /bin/bash does not get the path right to ABYSS currently so use /bin/sh
   say SCRIPT '#!/bin/sh';

   say SCRIPT "#SBATCH -J ${name}${kmer}";
   say SCRIPT "#SBATCH --mem=${mem_per_core}G";
   say SCRIPT "##SBATCH --partition=ksu-gen-highmem.q";
   say SCRIPT "#SBATCH --time=${maxtime}-23:00:00";
   say SCRIPT "#\$ -e ${kmer_dir}/kmer${kmer}.error";
   say SCRIPT "#SBATCH --output=${kmer_dir}/kmer${kmer}.out";
   say SCRIPT "#\$ -M $email";      # To email you when done
   say SCRIPT '#$ -m abe'; #to email when statuses change
#  say SCRIPT "#\$ ${exclusive}";
   say SCRIPT "#\$ ${priority}";
#   say SCRIPT "#\$ -q \'${queue}\'";

   if( $nodes == 1 ) {    # Single computer multi-core run
      say SCRIPT "#SBATCH --ntasks-per-node=${cores}\n";
   } else {
      say SCRIPT "#\$ -pe mpi-fill ${cores}\n";
   }


   say "Abyss executible is ".$abyss_exec;

   say SCRIPT "cd $cwd";
   say SCRIPT 'echo "SLURM_JOB_ID = $SLURM_JOB_ID"';
   say SCRIPT 'echo "SLURM_NPROCS = $SLURM_NPROCS"';
   say SCRIPT 'echo "np = ${np}"';
   say SCRIPT "hostname\n";

   say SCRIPT "module load ABySS\n";

	#Write Abyss driver script
	say SCRIPT "abyss-pe name=${name}${kmer} k=${kmer} " . 'np=$SLURM_NPROCS ' .
 		  "$single_names $paired_names $mp_names $single_files $paired_files $mp_files -C $kmer_dir";
 		  #WTF IS THIS STUFF
 		  #.
 		  #"${name}${kmer}-1.fa";



   close( SCRIPT );


#qsub if we are on Beocat or bail otherwise
   if ( `hostname` eq "eos" || "selene")
   {
#   		qsub will be instructed to wait for file chunking to finish
      #`qsub -hold_jid chunk1,chunk2 ${scriptfile}`;
    	`sbatch ${scriptfile}`;
      say "Job submitted\n";
   } else {
      say "Bailing since we are not on Beocat!";
      exit(0);
   }

   $kmer=$kmer+$kmer_increment;
}





sub load_data
{
	#Subroutine for testing loading/parsing of run data from a structured file
	#Requires that files are chunked by chunk_pair
	#Uses chunk_pair log style
	#File structure
	#[LIBRARY_NAME]
	#Type: (Single/Paired/Mate-pair) single is single read, paired is standard short insert, mate-pair is long span read pairs
	#Read1: (list of chunked files)
	#Read2: (list of chunked files)

	#Returns reference to hash holding assembly data

	my $file=shift; #File with reads in structured format
	my $dir=shift; #Adds directory to files for absolute path resolution

	open(my $fh, '<', $file);

	my %library_set; #holds data for a library
	my $library_name;
	while(<$fh>){
		chomp;

		unless (%library_set){
			#say "skipping first time";
		}

		#First look for library name
		if( /\[(.*)\]/ ){
			$library_name=$1;
		}
		elsif( /Type\: (.*)/ ){
			$library_set{$library_name}{'Type'}=$1;
		}
		elsif( /Read1\: (.*)/ ){
			my $r1=add_path($1, $cwd);
			$library_set{$library_name}{'Read1'}=$r1;
		}
		elsif( /Read2\: (.*)/ ){
			my $r2=add_path($1, $cwd);
			$library_set{$library_name}{'Read2'}=$r2;
		}
			elsif( /Singles\: (.*)/ ){
			my $s=add_path($1, $cwd);
			$library_set{$library_name}{'Singles'}=$s;
		}

	}
	return \%library_set;
}

sub add_path
{
	my $data=shift; #Get text that looks like 'Read1: file1 file2 ...' and add directory to it
	my $dir=shift;

	my @file=split(' ', $data);

	my @dir_files;
	foreach my $f (@file){
		$f=$dir.$f;
		#say $f;
		push(@dir_files, $f);
	}

	my $dir_data=join(' ', @dir_files);
	return $dir_data;
}




__END__

=head1 abyss_master3.pl

This script is a driver for running Abyss de novo (link) multi-k assemblies. A key limitation of Abyss is that reading/hashing files is single threaded and very inefficient if the number of files does not optimally match the number of CPUs being used for assembly.

Thus, before running this script, use chunk_pair4.pl to prepare your files into chunks that match the number of CPU cores on your machine. Abyss runs fastest when fastq files are <2GB compressed or <6GB uncompressed.

This script takes the output logs from chunk_pair4.pl to find and classify file types. After running chunk_pair4.pl on each library, create a file that logs all the files, their names and types. The format for the file is as follows:

[LIBRARY_NAME]
Type: (Single/Paired/Mate-pair) single is single read, paired is standard short insert, mate-pair is long span read pairs
Read1: (list of chunked files)
Read2: (list of chunked files)

For Single read libraries, omit the Read2 tag. LIBRARY_NAME is a unique identifier for the library type and must be a single word.

All files must be in the current directory and not have path information. The script determines its current location and add the path location information.

This script is designed for use either on local machines, or on a HPC compute cluster. At the moment it is configured for Beocat using the Sun Grid Engine. It could be easily adapted to other schedulers.

USAGE: abyss-master3.pl --data_file FILE

Options
--outdir	Specify output directory for assemblies, default abyss-assembly

--name		Specify the name for the assembly, default 'abyss'

--kstart	kmer start value, default 21

--kstop		kmer end value, default 31

--increment	kmer increment value, default 5

--nodes		Number of nodes to use, default 1

--cores		Number of cores to use, default 16

--maxtime	Maximum run time on the job scheduler in days, default 1

--mem		Total memory to use in GB. Note that for SGE, the memory per core is the total memory divided by the number of cores. Default 256 (16 per core)

--queue		**THIS OPTION HAS BEEN DISABLED** Specify a queue to use, default 'hero45'

--email		email address to notify you of SGE changes
