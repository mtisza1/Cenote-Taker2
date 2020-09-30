#!/bin/perl
use strict;
use warnings;
use File::Find;
use Cwd qw(cwd);

my$Drctry = cwd;
my@BlastOutput;


find(\&Wanted, $Drctry);

sub Wanted{
    if ($_ =~ m/.rotate.blastp.out/){
		push(@BlastOutput, $File::Find::name);
	}
}

foreach my$BlastPReport (@BlastOutput){
	print "$BlastPReport\n";
	my$Tbl = $BlastPReport;
	$Tbl =~ s/(\w+).rotate.blastp.out/$1.BLASTP.tbl/;
	my$Label = $1;
	open(MYFILE, "$BlastPReport") or die "$!";
	open(OUTPUT, ">$Tbl");
	print OUTPUT ">Feature $Label Table1";
	while(defined(my$line = <MYFILE>)){
		if($line =~ m/Query= ([^\s]+) \[(\d+) - (\d+)\]/){
			print OUTPUT "\n$2\t$3\tCDS";
			print OUTPUT "\n\t\t\tprotein_id\tlcl|$1";
		}
		elsif($line =~ m/\Q***** No hits found *****\E/){
			print OUTPUT "\n\t\t\tproduct\thypothetical protein";
			print OUTPUT "\n\t\t\tnote\tNo BLASTP hit";			
		}
		elsif($line =~ m/\>(.+\.\d) (.+) (\[|\.|\n)/){
			print OUTPUT "\n\t\t\tproduct\t$2";
			print OUTPUT "\n\t\t\tinference\tsimilar to AA sequence:INSD:$1";			

		}
	}
}