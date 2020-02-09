#!/bin/bash

git clone https://github.com/mtisza1/Cenote-Taker2.git

cd Cenote-Taker2
#chmod +x irf?
git clone --recursive https://github.com/deprekate/PHANOTATE.git
cd PHANOTATE; make
cd ..
wget http://last.cbrc.jp/last-1047.zip
unzip last-1047.zip
cd last-1047/
make
#### set apc script to call last correctly

conda env create --file cenote-taker2_env.yml

#need a way to run updatetaxonomy for krona
#/gpfs/gsfs7/users/tiszamj/python/env/Cenote-Taker2/opt/krona/updateTaxonomy.sh

conda activate Cenote-Taker2
# git hhsuite
#### i may want to "fork" PHANOTATE and hhsuite
#wget hmmer DBs
wget https://zenodo.org/record/3660539/files/hmmscan_DBs.tgz
tar -xvf hmmscan_DBs.tgz
rm hmmscan_DBs.tgz
#wget BLAST DBs
wget https://zenodo.org/record/3660538/files/blast_DBs.tgz
tar -xvf blast_DBs.tgz
rm blast_DBs.tgz
#wget hhsuite DBs
## NCBI_CD
wget https://zenodo.org/record/3660537/files/NCBI_CD_hhsuite.tgz
tar -xvf NCBI_CD_hhsuite.tgz
rm NCBI_CD_hhsuite.tgz
## pfam_32
mkdir pfam_32_db && cd pfam_32_db
wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_32.0.tar.gz
tar -xvf pfamA_32.0.tar.gz
rm pfamA_32.0.tar.gz
cd ..
## pdb70_latest
mkdir pdb70 && cd pdb70
wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz
tar -xvf pdb70_from_mmcif_latest.tar.gz
rm pdb70_from_mmcif_latest.tar.gz
cd ..
#wget cdd rpsblast db
mkdir cdd_rps_db && cd cdd_rps_db
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
tar -xvf Cdd_LE.tar.gz
rm Cdd_LE.tar.gz

#format DBs?

#python setup.py build; python setup.py install

#Cenote-Taker2 -h