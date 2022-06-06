#!/bin/perl

# AUTHOR: Joseph Fass
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2015 The Regents of The University of California, Davis Campus.
# All rights reserved.

use Getopt::Std;

my $usage = "\nusage: $0 [options] <fasta>\n\n".
            "Have you assembled [a] [p]erfect [c]ircle (e.g. plasmid, bacterial chromosome)? ".
            "This script tests each sequence in the supplied fasta file for self-overlap. ".
            "If overlap is found, the 5' copy of the overlapping sequence is trimmed, ".
            "the ends joined, and the join moved (\"permuted\") to the middle of the output sequence. ".
            "Join location is appended to the sequence header, to allow confirmation. ".
            "If a multi-fasta file is provided, each sequence is tested for circularity ".
            "and output separately. ".
            "If reads are provided (e.g. PacBio corrected.fastq), ".
            "BWA MEM is used to align reads to spliced sequence, for validation.\n".
            "The binaries \'lastdb\', \'lastal\', and (with -r) ".
            "\'bwa\' and \'samtools\' must be in your PATH.\n".
            "\n".
            "-b <text>\tbasename for output files (alphanumeric and underscore characters)\n".
            "-r <fastq>\thigh accuracy long reads to use for validation\n".
            "-c <text>\tlastdb command\n".
            "-d <text>\tlastal command\n".

            "\n";
our($opt_b,$opt_r,$opt_c,$opt_d);
getopts('r:b:c:d:') or die $usage;
if (!(defined($opt_b)) or !($opt_b=~m/^\w+$/)) {$opt_b="permuted"} 
$fname = $opt_b;  # works better in interpreted quotes
my $assemblyfile = shift or die $usage;

# pull in (single-line) assembled molecules
open ASM, "<$assemblyfile";
while ($line = <ASM>) {
    if ($line =~ m/^>/) {
	push @contigs, $block if !($block eq '');  # push 2-line contig onto list of contigs
	$block = $line."\n";  # adds another newline, because it'll be chomped below
    }
    else {
	chomp $block;  # remove newline before adding more sequence
	chomp $line;  # remove newline on current (sequence) line
	$block .= $line."\n";  # block always ends in a single newline
    }
     }
push @contigs, $block;  # push last sequence block onto list
close ASM;

# loop through each contig, perform LAST self-alignment, trim, output
$n = 0;  # keeps track of sequence count
while (@contigs) {
    #print "contig\t";     # debug
    #print $#contigs + 1;  # debug
    #print "\n";           # debug
    $block = shift(@contigs);
    $n++;
    open SEQ, ">temp.apc.fa";
    print SEQ $block;
    close SEQ;
    print "formatting LAST db ... ";
    $command = "$opt_c temp.apc temp.apc.fa";  # format lastdb for current contig
    system($command);
    print "running LAST self-alignment ... ";
    $command = "$opt_d -s 1 -x 300 -f 0 -T 1 temp.apc temp.apc.fa > temp.apc.lastal";  # self align current contig
    system($command);
    $command = "$opt_d -s 1 -x 300 -f BlastTab -T 1 temp.apc temp.apc.fa | grep -v \"^#\" > temp.apc.blasttab";  # self align current contig
    system($command);    
    print "done\n";
    # pull in output of current contig's LAST self-alignment
    open LAST, "<temp.apc.lastal";
    undef(@alignment);
    while ($line = <LAST>) {
	if (!($line =~ m/^#/)) {
	    push @alignment, $line;
	}
    }
    close LAST;
    # process alignment lines; should be three, with first being full self-alignment
    # second and third lines should be near identical end-overlaps
    if ($#alignment + 1 == 3) {
	# check that first is full self alignment
	#print "three alignments in LAST output!\n";  # debug
	# check that 2nd and 3rd are near identical overlaps
	@overlap5prime = split("\t", $alignment[1]);
	@overlap3prime = split("\t", $alignment[2]);
	#print $overlap5prime[0]."\n";  # debug
	#print $overlap3prime[0]."\n";  # debug
	if ($overlap5prime[0] == $overlap3prime[0]) {
		#$ONEQ = system("sed 's/ \\+/ /g' temp.apc.maf | grep \"^s \" | cut -d \" \" -f7 | head -n3 | tail -n1");
		#$TWOQ = system("sed 's/ \\+/ /g' temp.apc.maf | grep \"^s \" | cut -d \" \" -f7 | head -n4 | tail -n1");
		my $IDENT = `cut -f3 temp.apc.blasttab | head -n2 | tail -n1 `;
		#print $IDENT."\n";

		my $HUND=100.00;
		#print $HUND."\n";
		$EQQ = $IDENT eq $HUND;
		#print $EQQ."\n";
		if ($IDENT == $HUND) {
	    $start5prime = $overlap5prime[2] + 1;  # 0-based!
	    $end5prime = $start5prime + $overlap5prime[3] - 1;
	    $start3prime = $overlap3prime[2] + 1;  # 0-based!
	    $end3prime = $start3prime + $overlap3prime[3] - 1;
	    print "Overlap detected between  ";
	    print $overlap5prime[1].":".$start5prime."-".$end5prime;
	    print "  and  ";
	    print $overlap3prime[1].":".$start3prime."-".$end3prime."\n";

	    # trim one overlap, then permute halfway
	    $block =~ m/^(.*\n)(.*)$/;
	    $header = $1;
	    chomp $header;
	    $sequence = $2;
	    chomp $sequence;
	    $trimmedSeq = substr($sequence, $overlap5prime[3]);  # 0-, 1-based trickery
	    $mid = int( length($trimmedSeq) / 2 );  # midpoint, to permute around
	    $newSeq = substr($trimmedSeq, $mid) . substr($trimmedSeq, 0, $mid) . "\n";
	    open OUT, ">$fname.$n.fa";
	    chomp $block;
	    print OUT $header."|join_near_$mid\n";
	    print OUT $newSeq."\n";
	    close OUT;
	    open OUT2, ">$overlap5prime[1].DTR.tbl";
	    print OUT2 ">Feature $overlap5prime[1] Table1";
	    print OUT2 "\n".$start5prime."\t".$end5prime."\trepeat_region";
	    print OUT2 "\n\t\t\trpt_type\tDTR";
	    print OUT2 "\n".$start3prime."\t".$end3prime."\trepeat_region";
	    print OUT2 "\n\t\t\trpt_type\tDTR\n";
	    close OUT2;
	    if (defined($opt_r)) {
		# use bwa mem to align reads to permuted sequence
		$command = "bwa index $fname.$n.fa 2> $fname.$n.aln.log";
		print "indexing permuted sequence ... ";
		system($command);
		$command = "bwa mem -M $fname.$n.fa $opt_r 2>> $fname.$n.aln.log | gzip > $fname.$n.sam.gz";
		print "aligning reads to permuted sequence ... ";
		system($command);
		$command = "samtools view -uhS $fname.$n.sam.gz 2>> $fname.$n.aln.log | samtools sort - $fname.$n 2>> $fname.$n.aln.log";
		print "sorting BAM file of aligned reads ... ";
		system($command);
		$command = "samtools index $fname.$n.bam 2>> $fname.$n.aln.log";
		print "indexing BAM file ... ";
		system($command);
		print "done\n";
		}
		}

	}
	else {
	    system("mv temp.apc.lastal apc_aln_$n.txt");
	    print "irregular LAST results ... check local file apc_aln_$n.txt\n";
	}
    }
    else {
	system("mv temp.apc.lastal apc_aln_$n.txt");
	print "irregular LAST results ... check local file apc_aln_$n.txt\n";
    }
}

# clean up temp files!
system("rm temp.apc*");

