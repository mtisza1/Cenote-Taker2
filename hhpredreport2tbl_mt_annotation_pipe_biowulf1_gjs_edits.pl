#!/bin/perl
use strict;
use warnings;
use File::Find;
use Cwd qw(cwd);

my$Drctry = cwd;
my@HHPredOutput;


find(\&Wanted, $Drctry);

sub Wanted{
    if ($_ =~ m/.rotate.out_all.hhr/){
		push(@HHPredOutput, $File::Find::name);
	}
}

foreach my$HHPredReport (@HHPredOutput){
	print "$HHPredReport\n";
	my$Tbl = $HHPredReport;
	$Tbl =~ s/(\w+).rotate.out_all.hhr/$1.HH.tbl/;
	my$Label = $1;
	open(MYFILE, "$HHPredReport") or die "$!";
	open(OUTPUT, ">$Tbl");
	print OUTPUT ">Feature $Label Table1";
	my $some_bool = 0;
  my $last_record;
	while(defined(my$line = <MYFILE>)){
		if($line =~ m/Query         (\w+) \[(\d+) - (\d+)\]/){
      if ($some_bool == 0 && defined($last_record)) {
        print OUTPUT $last_record;
      }        
        print OUTPUT "\n$2\t$3\tCDS";
        print OUTPUT "\n\t\t\tprotein_id\tlcl|$1";

        undef $last_record;
        $some_bool = 0;
      
		}
		elsif ($line =~ /No \d+\n/){
          next if $some_bool == 1;
      		my @data = split(/\|/, $line);
      		my $nextline = <MYFILE>;
          if ($nextline =~ m/(\>\S+) (.*(hypothetical|uncharacterized|uncharacterised|expressed protein|genomic scaffold|uncultured bacterium|genome sequencing data|predicted protein|genome sequencing data).*)/i) {
            $last_record = "\n\t\t\tproduct\t$2\n\t\t\tinference\tsimilar to AA sequence:UniProtKB:$1";
          }
      		elsif($nextline =~ m/(\>\S+) (.+)/){
      			print OUTPUT "\n\t\t\tproduct\t$2";
            print OUTPUT "\n\t\t\tinference\tsimilar to AA sequence:UniProtKB:$1";
            $some_bool = 1;
      		}
		}
		elsif ($line =~ /^ No Hit/) {
      		my @data = split(/\|/, $line);
      		my $nextline = <MYFILE>;
      		if($nextline =~ m/^\s*$/){
			print OUTPUT "\n\t\t\tproduct\thypothetical protein";
      print OUTPUT "\n\t\t\tnote\tNo HHsearch hits";

      }
		}

	}
      if ($some_bool == 0 && defined($last_record)) {
        print OUTPUT $last_record;
      }
}
