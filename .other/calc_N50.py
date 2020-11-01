#!/usr/bin/perl

use warnings;
use strict;
use 5.010;

my $infile=$ARGV[0];
open(IN, "<$infile");

my $length=0; #Stores length of each sequence
my $total=0;
my @contig_lengths;

while(<IN>)
{
	if(/^>/ || eof) #If this is a sequence ID, we need to dump length and total to @contig_lengths
	{
		if($length>0) #If there is a previous length value from the previous contig
		{
			$total+=$length; #Add $length to the total
			push @contig_lengths, $length;
		}
		$length=0; #Reset length for the next contig
	}
	else
	{
		s/\s//g;
		$length+=length($_);
	}	
}

#Now calculate the N50

@contig_lengths=sort{$b<=>$a} @contig_lengths; #Sort the contig lengths

my $current_length;
my $count=0;
foreach(my $j=0; $count<$total/2; $j++)
{
		$count+=$contig_lengths[$j];
		$current_length=$contig_lengths[$j];
}
say $current_length;
	

