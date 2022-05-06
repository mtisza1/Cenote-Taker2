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
cenote_script_path = os.getcwd()   
print(cenote_script_path) 

parser = argparse.ArgumentParser(description='Update databases associated with Cenote-Taker2. HMM (hmmer) databases: updated June 16th, 2021. CDD (hhsuite) database: not updated from original. PFAM (hhsuite) database: not updated from original. PDB (hhsuite) database: updated periodically by PDB. CDD (rpsblast) database: not updated from original. Taxonomy BLAST (protein) databases: updated May 6th, 2022. Taxdump database is updated by NCBI periodically.')

optional_args = parser.add_argument_group('Use options to pick databases to update.')

optional_args.add_argument("--hmm", dest="HMM_DB", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--protein", dest="PROTEIN", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--rps", dest="RPS", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--taxdump", dest="TAXDUMP", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--hhCDD", dest="HHCDD", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--hhPFAM", dest="HHPFAM", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')
optional_args.add_argument("--hhPDB", dest="HHPDB", type=str2bool, default=False, help=' Default: False -- choose: True -or- False')

args = parser.parse_args()

if str(args.PROTEIN) == "True":
    print ("running taxonomy BLAST (protein) database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/blast_DBs/'])
    isExist = os.path.exists(str(cenote_script_path) + '/blast_DBs')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/blast_DBs/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'https://zenodo.org/record/6525617/files/blast_DBs.tar.gz'])
    subprocess.call(['tar', '-xvf', 'blast_DBs.tar.gz'])
    subprocess.call(['rm', '-f', 'blast_DBs.tar.gz'])

if str(args.HMM_DB) == "True":
    print ("running HMM database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/hmmscan_DBs/'])
    isExist = os.path.exists(str(cenote_script_path) + '/hmmscan_DBs')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/hmmscan_DBs/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'https://zenodo.org/record/4966268/files/hmmscan_DBs.tgz'])
    subprocess.call(['tar', '-xvf', 'hmmscan_DBs.tgz'])
    subprocess.call(['rm', '-f', 'hmmscan_DBs.tgz'])

if str(args.RPS) == "True":
    print ("running RPSBLAST CDD database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/cdd_rps_db/'])
    isExist = os.path.exists(str(cenote_script_path) + '/cdd_rps_db')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/cdd_rps_db/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz'])
    subprocess.call(['tar', '-xvf', 'Cdd_LE.tar.gz', '--directory', str(cenote_script_path) + '/cdd_rps_db'])
    subprocess.call(['rm', '-f', 'Cdd_LE.tar.gz'])

if str(args.TAXDUMP) == "True":
    print ("running TAXDUMP database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/taxdump/'])
    isExist = os.path.exists(str(cenote_script_path) + '/taxdump')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/taxdump/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'])
    subprocess.call(['tar', '-xvf', 'taxdump.tar.gz', '--directory', str(cenote_script_path) + '/taxdump'])
    subprocess.call(['rm', '-f', 'taxdump.tar.gz'])

if str(args.HHCDD) == "True":
    print ("running hhsuite CDD database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/NCBI_CD/'])
    isExist = os.path.exists(str(cenote_script_path) + '/NCBI_CD')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/NCBI_CD/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'https://zenodo.org/record/3660537/files/NCBI_CD_hhsuite.tgz'])
    subprocess.call(['tar', '-xvf', 'NCBI_CD_hhsuite.tgz'])
    subprocess.call(['rm', '-f', 'NCBI_CD_hhsuite.tgz'])

if str(args.HHPFAM) == "True":
    print ("running PFAM database update/install")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/pfam_32_db/'])
    isExist = os.path.exists(str(cenote_script_path) + '/pfam_32_db')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/pfam_32_db/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_32.0.tar.gz'])
    subprocess.call(['tar', '-xvf', 'pfamA_32.0.tar.gz', '--directory', str(cenote_script_path) + '/pfam_32_db'])
    subprocess.call(['rm', '-f', 'pfamA_32.0.tar.gz'])

if str(args.HHPDB) == "True":
    print ("running PDB database update/install. This could take around 2 hours.")
    subprocess.call(['rm', '-r', '-f', str(cenote_script_path) + '/pdb70/'])
    isExist = os.path.exists(str(cenote_script_path) + '/pdb70')
    if not isExist:
        os.makedirs(str(cenote_script_path) + '/pdb70/', exist_ok=True)
    subprocess.call(['wget', '--directory-prefix=' + str(cenote_script_path), 'http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz'])
    subprocess.call(['tar', '-xvf', 'pdb70_from_mmcif_latest.tar.gz', '--directory', str(cenote_script_path) + '/pdb70'])
    subprocess.call(['rm', '-f', 'pdb70_from_mmcif_latest.tar.gz'])
