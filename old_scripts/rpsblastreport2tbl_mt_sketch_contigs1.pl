#!/bin/perl
use strict;
use warnings;
use File::Find;
use Cwd qw(cwd);

my$Drctry = cwd;
my@BlastOutput;


find(\&Wanted, $Drctry);

sub Wanted{
    if ($_ =~ m/.AA.fasta.rpsblast1.out/){
		push(@BlastOutput, $File::Find::name);
	}
}

foreach my$RPSBlastReport (@BlastOutput){
	print "$RPSBlastReport\n";
	my$Tbl = $RPSBlastReport;
	$Tbl =~ s/(\w+).AA.fasta.rpsblast1.out/$1.NT.tbl/;
	my$Label = $1;
	open(MYFILE, "$RPSBlastReport") or die "$!";
	open(OUTPUT, ">$Tbl");
	print OUTPUT ">Feature $Label Table1";
	while(defined(my$line = <MYFILE>)){
		if($line =~ m/Query= (\w+) \[(\d+) - (\d+)\]/){
			print OUTPUT "\n$2\t$3\tCDS";
			print OUTPUT "\n\t\t\tprotein_id\tlcl|$1";

		}
		elsif($line =~ m/\Q***** No hits found *****\E/){
			print OUTPUT "\n\t\t\tproduct\thypothetical protein";
			print OUTPUT "\n\t\t\tnote\tNo RPS-BLAST hit";			
		}
		elsif($line =~ m/\>(\S+) (\Qcd\E\d+), (\S+), (.+)(\[|\.|\n)/){
			print OUTPUT "\n\t\t\tproduct\t$3";
			print OUTPUT "\n\t\t\tinference\tprotein motif:$1";			
		}		
		elsif($line =~ m/\>(\S+) (\S+), (\S+), (.+)(\[|\.|\n)/){
			print OUTPUT "\n\t\t\tproduct\t$4";
			print OUTPUT "\n\t\t\tinference\tprotein motif:$1";			
		}
	}
}