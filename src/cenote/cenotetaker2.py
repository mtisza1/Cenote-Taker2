#!/usr/bin/env python

import argparse
import sys, os
import subprocess
from pathlib import Path

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
pathname = os.path.dirname(__file__)  
cenote_script_path = os.path.abspath(pathname)      
print("this script dir: " + str(cenote_script_path))

parentpath = Path(pathname).parents[1]

__version__ = "2.2.0"

#def runner():
parser = argparse.ArgumentParser(description='Cenote-Taker 2 is a pipeline for virus discovery and thorough annotation of viral contigs and genomes. Visit https://github.com/mtisza1/Cenote-Taker2#use-case-suggestionssettings for suggestions about how to run different data types and https://github.com/mtisza1/Cenote-Taker2/wiki to read more. Version ' + str(__version__))

required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for Cenote-Taker2 ')

required_args.add_argument("-c", "--contigs", dest="original_contigs", type=str, required=True, help='Contig file with .fasta extension in fasta format. Each header must be unique before the first space character')
required_args.add_argument("-r", "--run_title", dest="run_title", type=str, required=True, help='Name of this run. A directory of this name will be created. Must be unique from older runs or older run will be renamed. Must be less than 18 characters, using ONLY letters, numbers and underscores (_)')

required_args.add_argument("-p", "--prune_prophage", dest="PROPHAGE", type=str2bool, required=True, help='True or False. Attempt to identify and remove flanking chromosomal regions from non-circular contigs with viral hallmarks (True is highly recommended for sequenced material not enriched for viruses. Virus-enriched samples probably should be False (you might check enrichment with ViromeQC). Also, please use False if --lin_minimum_hallmark_genes is set to 0)')
required_args.add_argument("-m", "--mem", dest="MEM", type=int, required=True, help='example: 56 -- Gigabytes of memory available for Cenote-Taker2. Typically, 16 to 32 should be used. Aim for at least 1/2 the value of \'-\' ')
required_args.add_argument("-t", "--cpu", dest="CPU", type=int, required=True, help='Example: 32 -- Number of CPUs available for Cenote-Taker2. Approximately 32 CPUs should be used moderately sized metagenomic assemblies. For large datasets, increased performance can be seen up to 120 CPUs. Fewer than 16 CPUs will work in theory, but could extend the length of the run. See GitHub repo for suggestions. ')




optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for Cenote-Taker2. Most of which are important to consider!!! GenBank typically only accepts genome submission with ample metadata. See https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#ModifiersPage for more information on GenBank metadata fields')

optional_args.add_argument('--version', action='version', version=str(__version__))
optional_args.add_argument("-am", "--annotation_mode", dest="ANNOTATION_MODE", type=str2bool, default="False", help='Default: False -- Annotate sequences only (skip discovery). Only use if you believe each provided sequence is viral')
optional_args.add_argument("--template_file", dest="template_file", type=str, default=str(parentpath) + '/dummy_template.sbt', help='Template file with some metadata. Real one required for GenBank submission. Takes a couple minutes to generate: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ ')
optional_args.add_argument("--reads1", dest="F_READS", type=os.path.abspath, default='no_reads', help=' Default: no_reads -- ILLUMINA READS ONLY: First Read file in paired read set - OR - read file in unpaired read set - OR - read file of interleaved reads. Used for coverage depth determination.')
optional_args.add_argument("--reads2", dest="R_READS", type=os.path.abspath, default='no_reads', help='Default: no_reads -- ILLUMINA READS ONLY: Second Read file in paired read set. Disregard if not using paired reads. Used for coverage depth determination.')
optional_args.add_argument("--minimum_length_circular", dest="circ_length_cutoff", type=int, default='1000', help='Default: 1000 -- Minimum length of contigs to be checked for circularity. Bare minimun is 1000 nts')
optional_args.add_argument("--minimum_length_linear", dest="linear_length_cutoff", type=int, default='1000', help='Default: 1000 -- Minimum length of non-circualr contigs to be checked for viral hallmark genes.')
optional_args.add_argument("-db", "--virus_domain_db", dest="virus_domain_db", type=str, default='virion', help='default: virion -- \'standard\' database: all virus (DNA and RNA) hallmark genes (i.e. genes with known function as virion structural, packaging, replication, or maturation proteins specifically encoded by virus genomes) with low false discovery rate. \'virion\' database: subset of \'standard\', hallmark genes encoding virion structural proteins, packaging proteins, or capsid maturation proteins (DNA and RNA genomes) with LOWEST false discovery rate. \'rna_virus\' database: For RNA virus hallmarks only. Includes RdRp and capsid genes of RNA viruses. Low false discovery rate.')
optional_args.add_argument("--lin_minimum_hallmark_genes", dest="LIN_MINIMUM_DOMAINS", type=int, default='1', help='Default: 1 -- Number of detected viral hallmark genes on a non-circular contig to be considered viral and recieve full annotation. WARNING: Only choose \'0\' if you have prefiltered the contig file to only contain putative viral contigs (using another method such as VirSorter or DeepVirFinder), or you are very confident you have physically enriched for virus particles very well (you might check with ViromeQC). Otherwise, the duration of the run will be extended many many times over, largely annotating non-viral contigs, which is not what Cenote-Taker2 is meant for. For unenriched samples, \'2\' might be more suitable, yielding a false positive rate near 0. ')
optional_args.add_argument("--circ_minimum_hallmark_genes", dest="CIRC_MINIMUM_DOMAINS", type=int, default='1', help='Default:1 -- Number of detected viral hallmark genes on a circular contig to be considered viral and recieve full annotation. For samples physically enriched for virus particles, \'0\' can be used, but please treat circular contigs without known viral domains cautiously. For unenriched samples, \'1\' might be more suitable. ')
optional_args.add_argument("--known_strains", dest="handle_knowns", type=str, default='do_not_check_knowns', help='Default: do_not_check_knowns -- do not check if putatively viral contigs are highly related to known sequences (via MEGABLAST). \'blast_knowns\': REQUIRES \'--blastn_db\' option to function correctly. ')
optional_args.add_argument("--blastn_db", dest="BLASTN_DB", type=str, default='none', help='Default: none -- Set a database if using \'--known_strains\' option. Specify BLAST-formatted nucleotide datase. Probably, use only GenBank \'nt\' database, \'nt viral\', or a subset therof, downloaded from ftp://ftp.ncbi.nlm.nih.gov/ Headers must be GenBank record format')
optional_args.add_argument("--enforce_start_codon", dest="ENFORCE_START_CODON", type=str2bool, default=False, help='Default: False -- For final genome maps, require ORFs to be initiated by a typical start codon? GenBank submissions containing ORFs without start codons can be rejected. However, if True,  important but incomplete genes could be culled from the final output. This is relevant mainly to contigs of incomplete genomes ')
optional_args.add_argument("-hh", "--hhsuite_tool", dest="HHSUITE_TOOL", type=str, default='hhblits', help=' default: hhblits -- hhblits will query PDB, pfam, and CDD to annotate ORFs escaping identification via upstream methods. \'hhsearch\': hhsearch, a more sensitive tool, will query PDB, pfam, and CDD to annotate ORFs escaping identification via upstream methods. (WARNING: hhsearch takes much, much longer than hhblits and can extend the duration of the run many times over. Do not use on large input contig files). \'no_hhsuite_tool\': forgoes annotation of ORFs with hhsuite. Fastest way to complete a run. ')
optional_args.add_argument("--crispr_file", dest="CRISPR_FILE", type=str, default='none', help='Tab-separated file with CRISPR hits in the following format: CONTIG_NAME HOST_NAME NUMBER_OF_MATCHES. You could use this tool: https://github.com/edzuf/CrisprOpenDB. Then reformat for Cenote-Taker 2')
optional_args.add_argument("--isolation_source", dest="isolation_source", type=str, default='unknown', help='Default: unknown -- Describes the local geographical source of the organism from which the sequence was derived')
optional_args.add_argument("--Environmental_sample", dest="Environmental_sample", type=str2bool, default=False, help='Default: False -- True or False, Identifies sequence derived by direct molecular isolation from an unidentified organism')
optional_args.add_argument("--collection_date", dest="collection_date", type=str, default='unknown', help='Default: unknown -- Date of collection. this format: 01-Jan-2019, i.e. DD-Mmm-YYYY')
optional_args.add_argument("--metagenome_type", dest="metagenome_type", type=str, default='unknown', help='Default: unknown -- a.k.a. metagenome_source')
optional_args.add_argument("--srr_number", dest="srr_number", type=str, default='unknown', help='Default: unknown -- For read data on SRA, run number, usually beginning with \'SRR\' or \'ERR\' ')
optional_args.add_argument("--srx_number", dest="srx_number", type=str, default='unknown', help='Default: unknown -- For read data on SRA, experiment number, usually beginning with \'SRX\' or \'ERX\' ')
optional_args.add_argument("--biosample", dest="biosample", type=str, default='unknown', help='Default: unknown -- For read data on SRA, sample number, usually beginning with \'SAMN\' or \'SAMEA\' or \'SRS\' ')
optional_args.add_argument("--bioproject", dest="bioproject", type=str, default='unknown', help='Default: unknown -- For read data on SRA, project number, usually beginning with \'PRJNA\' or \'PRJEB\' ')
optional_args.add_argument("--assembler", dest="ASSEMBLER", type=str, default='unknown_assembler', help='Default: unknown_assembler -- Assembler used to generate contigs, if applicable. Specify version of assembler software, if possible. ')
optional_args.add_argument("--molecule_type", dest="MOLECULE_TYPE", type=str, default='DNA', help='Default: DNA -- viable options are DNA - OR - RNA ')
optional_args.add_argument("--data_source", dest="DATA_SOURCE", type=str, default='original', help=' default: original -- original data is not taken from other researchers\' public or private database. \'tpa_assembly\': data is taken from other researchers\' public or private database. Please be sure to specify SRA metadata.  ')
optional_args.add_argument("--filter_out_plasmids", dest="FILTER_PLASMIDS", type=str2bool, default=True, help='Default: True -- True - OR - False. If True, hallmark genes of plasmids will not count toward the minimum hallmark gene parameters. If False, hallmark genes of plasmids will count. Plasmid hallmark gene set is not necessarily comprehensive at this time. ')
optional_args.add_argument("--scratch_directory", dest="SCRATCH_DIR", type=str, default="none", help='Default: none -- When running many instances of Cenote-Taker2, it seems to run more quickly if you copy the hhsuite databases to a scratch space temporarily. Use this argument to set a scratch directory that the databases will be copied to (at least 100GB of scratch space are required for copying the databases)')
optional_args.add_argument("--blastp", dest="BLASTP", type=str, default="no_blastp", help='Do not use this argument as of now. ')
optional_args.add_argument("--orf-within-orf", dest="ORF_WITHIN", type=str2bool, default="False", help='Default: False -- Remove called ORFs without HMMSCAN or RPS-BLAST hits that begin and end within other ORFs? True or False')
optional_args.add_argument("--cenote-dbs", dest="CENOTE_DBS", type=str, default=parentpath, help='Default: ' + str(parentpath) + ' -- If you downloaded and setup the databases in a non-standard location, specify path')
optional_args.add_argument("--wrap", dest="WRAP", type=str2bool, default="True", help='Default: True -- Wrap/rotate DTR/circular contigs so the start codon of an ORF is the first nucleotide in the contig/genome')
optional_args.add_argument("--hallmark_taxonomy", dest="HALLMARK_TAX", type=str2bool, default="False", help='Default: False -- Get hierarchical taxonomy information for all hallmark genes? This report (.hallmarks.taxonomy.out) is not considered in the final taxonomy call.')
optional_args.add_argument("--exact_dtrs", dest="EXACT", type=str2bool, default="False", help='Default: False -- Only record DTRs if each ends matches with 100 percent ID.')



args = parser.parse_args()

def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None

if not is_tool("samtools") and str(args.F_READS) != "no_reads" :
    print ("samtools is not found. Exiting.")
    quit()    
if not is_tool("minimap2") :
    print ("minimap2 is not found. Exiting.")
    quit()
#    if not is_tool("bioawk") :
#        print ("bioawk is not found. Exiting.")
#        quit()
if not is_tool("tRNAscan-SE") :
    print ("tRNAscan-SE is not found. Exiting.")
    quit()
if not is_tool("tbl2asn") :
    print ("table2asn is not found. Exiting.")
    quit()
if not is_tool("seqkit") :
    print ("seqkit is not found. Exiting.")
    quit()
#    if not is_tool("hhblits") :
#        print ("hhblits is not found. Exiting.")
#        quit()
if not is_tool("bedtools") :
    print ("bedtools is not found. Exiting.")
    quit()
if not is_tool("phanotate.py") :
    print ("phanotate is not found. Exiting.")
    quit()
if not is_tool("prodigal") :
    print ("prodigal is not found. Exiting.")
    quit()

reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]


if 'pyhmmer' not in installed_packages:
    print("pyhmmer not found in installed python packages. Exiting.")
    quit()

if 'numpy' not in installed_packages:
    print("numpy not found in installed python packages. Exiting.")
    quit()

if 'pandas' not in installed_packages:
    print("pandas not found in installed python packages. Exiting.")
    quit()

if 'biopython' not in installed_packages:
    print("biopython not found in installed python packages. Exiting.")
    quit()



subprocess.call(['bash', str(cenote_script_path) + '/new_cenote1.sh', 
                    str(args.original_contigs), str(args.run_title)])
