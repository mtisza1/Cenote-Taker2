#!/usr/bin/env python

import argparse
import sys, os
import subprocess

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
pathname = os.path.dirname(sys.argv[0])  
cenote_script_path = os.path.abspath(pathname)      
print(cenote_script_path) 

parser = argparse.ArgumentParser(description='unlimited_breadsticks is a pipeline for virus discovery and cursory annotation of viral contigs and genomes. Visit https://github.com/mtisza1/Cenote_Unlimited_Breadsticks for further description')

required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for unlimited_breadsticks ')

required_args.add_argument("-c", "--contigs", dest="original_contigs", type=str, required=True, help='Contig file with .fasta extension in fasta format - OR - assembly graph with .fastg extension. Each header must be unique before the first space character')
required_args.add_argument("-r", "--run_title", dest="run_title", type=str, required=True, help='Name of this run. A directory of this name will be created. Must be unique from older runs or older run will be renamed. Must be less than 18 characters, using ONLY letters, numbers and underscores (_)')
required_args.add_argument("-p", "--prune_prophage", dest="PROPHAGE", type=str2bool, required=True, help='True or False. -- Attempt to identify and remove flanking chromosomal regions from non-circular contigs with viral hallmarks (True is highly recommended for sequenced material not enriched for viruses. Virus enriched samples probably should be False (you might check actaul enrichment with ViromeQC). Also, please use False if --lin_minimum_hallmark_genes is set to 0)')
required_args.add_argument("-m", "--mem", dest="MEM", type=int, required=True, help='example: 56 --	Gigabytes of memory available for unlimited_breadsticks. Typically, 16 to 32 should be used. Lower memory will work in theory, but could extend the length of the run ')
required_args.add_argument("-t", "--cpu", dest="CPU", type=int, required=True, help='Example: 32 --	Number of CPUs available for unlimited_breadsticks. Typically, 32 CPUs should be used. For large datasets, increased performance can be seen up to 120 CPUs. Fewer than 16 CPUs will work in theory, but could extend the length of the run ')




optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for unlimited_breadsticks. Most of which are important to consider!!! ')


optional_args.add_argument("--minimum_length_circular", dest="circ_length_cutoff", type=int, default='1000', help='Default: 1000 -- Minimum length of contigs to be checked for circularity. Absolute minimun is 1000 nts')
optional_args.add_argument("--minimum_length_linear", dest="linear_length_cutoff", type=int, default='1000', help='Default: 1000 -- Minimum length of non-circualr contigs to be checked for viral hallmark genes.')
optional_args.add_argument("-db", "--virus_domain_db", dest="virus_domain_db", type=str, default='virion', help='default: virion -- \'standard\' database: all virus (DNA and RNA) hallmark genes (i.e. genes with known function as virion structural, packaging, replication, or maturation proteins specifically encoded by virus genomes) with very low false discovery rate. \'virion\' database: subset of \'standard\', hallmark genes encoding virion structural proteins, packaging proteins, or capsid maturation proteins (DNA and RNA genomes). \'rna_virus\' database: For RNA virus hallmarks only. Includes RdRp and capsid genes of RNA viruses. Low false discovery rate due to structural similarity between RdRp genes and e.g. transposon-encoded RT genes')
optional_args.add_argument("--lin_minimum_hallmark_genes", dest="LIN_MINIMUM_DOMAINS", type=int, default='1', help='Default: 1 -- Number of detected viral hallmark genes on a non-circular contig to be considered viral. ')
optional_args.add_argument("--circ_minimum_hallmark_genes", dest="CIRC_MINIMUM_DOMAINS", type=int, default='1', help='Default:1 -- Number of detected viral hallmark genes on a circular contig to be considered viral. ')
optional_args.add_argument("--filter_out_plasmids", dest="FILTER_PLASMIDS", type=str2bool, default=True, help='Default: True -- True - OR - False. If True, hallmark genes of plasmids will not count toward the minimum hallmark gene parameters. If False, hallmark genes of plasmids will count. Plasmid hallmark gene set is not necessarily comprehensive at this time. ')



args = parser.parse_args()

def is_tool(name):
	"""Check whether `name` is on PATH."""
	from distutils.spawn import find_executable
	return find_executable(name) is not None

if is_tool("prodigal") :
	print ("prodigal found")
else:
	print ("prodigal is not found. Exiting.")
	quit()

if is_tool("rpsblast") :
	print ("rpsblast found")
else:
	print ("rpsblast is not found. Exiting.")
	quit()
if is_tool("bioawk") :
	print ("bioawk found")
else:
	print ("bioawk is not found. Exiting.")
	quit()
if is_tool("hmmscan") :
	print ("hmmscan found")
else:
	print ("hmmscan is not found. Exiting.")
	quit()



subprocess.call(['bash', str(cenote_script_path) + '/unlimited_breadsticks_0.1.sh', str(args.original_contigs), str(args.run_title), str(args.circ_length_cutoff), str(args.linear_length_cutoff), str(args.virus_domain_db), str(args.LIN_MINIMUM_DOMAINS), str(args.PROPHAGE), str(args.FILTER_PLASMIDS), str(cenote_script_path), str(args.CIRC_MINIMUM_DOMAINS), str(args.MEM), str(args.CPU)])

