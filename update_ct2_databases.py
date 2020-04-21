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

parser = argparse.ArgumentParser(description='Update databases associated with Cenote-Taker2. HMM (hmmer) databases: updated April 21st, 2020. CDD (hhsuite) database: not updated from original. PFAM (hhsuite) database: not updated from original. PDB (hhsuite) database: not updated from original. CDD (rpsblast) database: not updated from original. Taxonomy (BLAST) databases: not updated from original.')

optional_args = parser.add_argument_group('Use options to pick databases to update.')

optional_args.add_argument("--hmm", dest="HMM_DB", type=str2bool, default=False, help=' Default: false... Set to: True -or- False')

if str(args.HMM_DB) == True:
	subprocess.call(['rm -r ', str(cenote_script_path) + '/hmmscan_DBs/'])
	subprocess.call(['wget --directory-prefix=', str(cenote_script_path) + ' https://zenodo.org/record/3759823/files/hmmscan_DBs.tgz'])
	subprocess.call(['tar -xvf hmmscan_DBs.tgz'])
	subprocess.call(['rm hmmscan_DBs.tgz'])

