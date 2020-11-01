#!/usr/bin/perl

#SBATCH --time=20:00:00

use strict;
use warnings;
use 5.010;
use Getopt::Long;
#use YAML;

#Revised version of Dave Turner's Chunk Pair to interface better with abyss_master2.pl
#Version 4 pre-computes the chunked file name and writes them to a log for use by abyss_master2.pl
#This allows chunking to continue on a job scheduler, the scripts written without having to wait for chunking to complete

# USAGE: chunk_pair4.pl --read1 read1.fastq.gz --read2 read2.fastq.gz --nfiles n --library_name name --library_type 

# This takes about 10 minutes per compressed GB on a slow machine.

# pair is 1 or 2 and sequences will be marked /1 or /2
# This will gunzip the file if necessary, add pairing, chunk it,
# and gzip it at the end.
# Returns array reference to file chunk list



#my $t_start;
#my $time;

# Process the command line arguments
my ($read1, $read2, $singles, $nfiles, $pairing, $logfile, $library_name, $library_type);
GetOptions(	"read1=s"	=>	\$read1,
			"read2=s"	=>	\$read2,
			"singles=s"	=>	\$singles,
			"nfiles=i"	=>	\$nfiles,
			"pairing=i"	=>	\$pairing,
			"logfile=s" =>	\$logfile,
			"library_name=s"	=>	\$library_name,
			"library_type=s"	=>	\$library_type) or die ("Error in command line arguments\n");

my %log; #hash to log chunking data
#Start logging process
$log{$library_name}{'Type'}=$library_type;


#Call subroutine to chunk files
if($read1){
	my @files_read1=chunk_file($read1, $nfiles, '1');
}
if($read2){
	my @files_read2=chunk_file($read2, $nfiles, '2');
}			
if($singles){
	my @files_singles=chunk_file($singles, $nfiles, undef);
}

#Write log file that can be used for abyss_master3.pl
#print Dump(\%log);
write_log(\%log, $logfile);			
			
			
			
			
			
#Sub takes three arguments, file in fastq format, number of files to split it into, and if it is paired.
#Returns list of files created				
sub chunk_file
{
	my $infile=shift; #Fastq file to process
	my $nfiles=shift; #Number of files to split into
	my $pairing=shift; #Pairing 1, 2 of undef
	my $log=shift; #Array ref to log array
	my @filenames; #Array of filenames to return


	if( ! -r $infile )
	{
	   say "Could not find file $infile";
	   exit(-1);
	}

	#This determines if pairing is being used and sets opposite pair number
	my $oppopair;
	my $pair_tag; #Tag for hash of read types
	if( not defined($pairing) ) { $pairing=0; $oppopair=0; $pair_tag='Singles' } #Set to zero if single reads
	if( $pairing == 1 ) { $oppopair = 2; $pair_tag='Read1' }
	if( $pairing == 2 ) { $oppopair = 1; $pair_tag='Read2' }
	#say "Pairing set to $pairing";



	my $nlines;
	my $outfilebase=$infile; #Outputfile base
	#Testing if gzipped or not
	if($outfilebase =~ /(\.fastq\.gz)/)
	{
		#say "gzipped";
		$nlines = `gunzip $infile --stdout | wc -l`;
		open( LINES, "-|", "gunzip $infile --stdout" )
		  or die "cannot open $infile for input";
	}
	elsif($outfilebase =~ /(\.fastq)/)
	{
		#say "not gzipped";
		$nlines = `cat $infile | wc -l`;
		open( LINES, "<", $infile )
		  or die "cannot open $infile for input";

	}	
	chomp($nlines);
	my $nseq = $nlines / 4;

	#$time = time() - $t_start;

	#printf(" %d sequences found in file %s took %7.3lf seconds\n",
	#	   $nseq, $infile, $time);

	# Round nseq down leaving extra for the last file
	$nseq -= $nseq % $nfiles; 
	my $nseq_per = $nseq / $nfiles;

	my $nlines_per = $nseq_per * 4;

	#say "$nlines lines for $nfiles files is $nlines_per lines/file";
	#say "The extra sequences will go in the last file\n";

	#Build filenames, put into @filenames array, write logfile of names
	for(my $i=1; $i<=$nfiles; $i++){
		push(@filenames, "${outfilebase}-${i}.gz")
	}
	#write_log(@filenames, $pairing, $logfile); #subroutine writes log file before processing it
	$log{$library_name}{$pair_tag}="@filenames"; #Log names of files in hash

	# Open the input file and initialize line and file_num counters
	my $n = 0;            # line counter
	my $file_num = 0;
	my $LINES;

	# Process each line by dumping to the open file, changing to the next as needed

	#$t_start = time();
	my $line;
	while( $line = <LINES> )
	{
	   if( $n == 0 )  # open a new output file
	   {
		  if ( $file_num < $nfiles ){
			if( $file_num > 0 ) { close( OUTFILE ); }
			$file_num ++;
		
			my $current_file=shift(@filenames);

			open( OUTFILE, "|-", "gzip > $current_file")
				or die "cannot open $current_file for output";
			#push(@filenames, "${outfilebase}-${file_num}.gz");
			$n = $nlines_per;
			say "Opening $current_file for output";
		  } 
		  else {
			 $n = $nlines - ( $nlines_per * $nfiles );  # These are the extra lines to put in the last file
		  }
	   }

	# It is a sequence header so handle pairing
	# Check to avoid redundant pair labeling, exit if opposite pair label is there.

	   if( ( $line =~ /^@/ ) && ( $n % 4 == 0 ) && ( $pairing eq '1' || '2' || '0') ) {
		  if( $line =~ /\/$oppopair$/ ) {
			 say "ERROR:  A line already has the opposite pairing";
			 exit(-1);
		  }		  
		  if( $line =~ /\/$pairing$/ ) {  # Already marked with pairing
			 print OUTFILE $line;
		  } 
		  elsif($pairing == 0){
			#If these are single end reads do not add pairing
			chomp($line);
			#$line=~s/\s/\_/g; #BJSCO Must remove spaces from sequence header
								#BJSCO OR remove the " 1:N:0:7" at the end of the line
								#BJSCO note that : separated elements are variable
			#The " 1:N:0:7" data at the end of the read is messing up paring. Remove it
			$line=~s/\s.*//;
			 print OUTFILE "${line}\n";
		  } 
		  else {
			 chomp($line);
			 #$line=~s/\s/\_/g; #BJSCO Must remove spaces from sequence header
								#BJSCO OR remove the " 1:N:0:7" at the end of the line
								#BJSCO note that : separated elements are variable
			#The " 1:N:0:7" data at the end of the read is messing up paring. Remove it
			$line=~s/\s.*//;
						
			 print OUTFILE "${line}/${pairing}\n";
			 #print OUTFILE "${line}\n";
		  }
	   } 
	   else {
		  print OUTFILE $line;
	   }
	   $n--;
	}

	close($infile);
	close OUTFILE; #Required otherwise last file may not close


}

#Sub writes a log file in the format abyss_master3 can use to run
#Takes a hash as an argument
sub write_log
{
	my $log=shift; #Reference to hash with log data
	my $log_file=shift; #File to write log to
	
	open(my $log_fh, '>', $logfile) or die "Log file not specified! $!";
	
	
	for my $library_name (sort {$log{$a} <=> $log{$b}} keys %{ $log }){
		say $log_fh '['.$library_name.']';
		say $log_fh 'Type: '.$log{$library_name}{'Type'};
		delete $log{$library_name}{'Type'}; #A bit hacky but need to print out "Type: " first, then reads
		for my $read_tag (sort keys %{ $log{$library_name} }) {
			say $log_fh "$read_tag: ".$log{$library_name}{$read_tag};
		}
	}
	close $log_fh;

}



__END__

=head1 chunk_pair5.pl

This script takes Illumina fastq files and splits them up into chunks to prepare them for assembly with abyss-master3.pl. This script writes a log of the files that are broken up into chunks that is compatible as input for abyss-master3.pl.

USAGE: 
For paired reads, or mate-pair reads
chunk_pair5.pl --read1 FILE1.fastq --read2 FILE2.fastq --nfiles n --library_name LIBRARY_NAME --library_type (Single|Paired|Mate-paired) --logfile LOGFILE.txt

OPTIONS:

--read1	Read1 fastq file

--read2	Read2 fastq file

--nfiles	number of file chunks to create

--logfile	Specify a different name for the logfile

--library_name	Unique name for the library 

--library_type	Library type, must be Single, Paired or Mate-paired 


