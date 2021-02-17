#!/bin/perl
use strict;
use warnings;
use File::Find;
use Cwd qw(cwd);

my$Drctry = cwd;
my@BlastOutput;


find(\&Wanted, $Drctry);

sub Wanted{
    if ($_ =~ m/.AA.rpsblast.out/){
		push(@BlastOutput, $File::Find::name);
	}
}

foreach my$RPSBlastReport (@BlastOutput){
	print "$RPSBlastReport\n";
	my$Tbl = $RPSBlastReport;
	$Tbl =~ s/(\w+).AA.rpsblast.out/$1.RPS_TABLE.txt/;
	my$Label = $1;
	open(MYFILE, "$RPSBlastReport") or die "$!";
	open(OUTPUT, ">$Tbl");
	while(defined(my$line = <MYFILE>)){
		if($line =~ m/Query= ([^\s]+) \[(\d+) - (\d+)\]/){
			print OUTPUT "$1\t$2\t$3\t";

		}
		elsif($line =~ m/\Q***** No hits found *****\E/){
			print OUTPUT "HYPO\thypothetical protein\n";
		}
		elsif($line =~ m/\>(\S+) (\Qcd\E\d+), (\S+), (.+)(\[|\.|\n)/){
			print OUTPUT "$2\t$4\n";
		}		
		elsif($line =~ m/\>(\S+) (\S+), (\S+), (.+)(\[|\.|\n)/){
			print OUTPUT "$2\t$4\n";
		}
	}
}