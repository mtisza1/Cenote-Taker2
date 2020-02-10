#!/bin/bash

git clone https://github.com/mtisza1/Cenote-Taker2.git

cd Cenote-Taker2
#chmod +x irf307.linux.exe
# cloning PHANOTATE
git clone --recursive https://github.com/deprekate/PHANOTATE.git
cd PHANOTATE; # git checkout version from feb 8 2020
make
cd ..


wget http://last.cbrc.jp/last-1047.zip
unzip last-1047.zip
cd last-1047/
make
cd ..

CONDA_BASE=$( conda info --base )
eval "$(conda shell.bash hook)"
conda env create --file cenote-taker2_env.yml
#conda create -n Cenote-Taker2 -c defaults -c bioconda -c AgBiome --no-channel-priority python=3.6 prodigal=2.6.3 BWA=0.7.17 samtools=1.3 mummer=3.23 circlator=1.5.5 blast=2.9.0 bioawk=1.0 entrez-direct=13.3 krona=2.7.1 hmmer=3.3 bowtie2=2.3.5 trnascan-se=2.0.5 bbtools tbl2asn=25.7 emboss=6.6.0 cmake numpy pandas matplotlib 
conda activate Cenote-Taker2

# getting hh-suite from github. The anaconda package of hh-suite causes conflicts
git clone https://github.com/soedinglab/hh-suite.git
cd hh-suite # git checkout version from feb 8 2020
mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
cd ..
cd ..


conda info --envs | sed 's/ \+/ /g' | if grep -q "Cenote-Taker2 \*" ; then 
	echo "Cenote-Taker2 loaded" ; 
else 
	echo "Cenote-Taker2 not loaded correctly" ;
	exit 
fi

KRONA_DIRE=$( which python | sed 's/bin\/python/opt\/krona/g' )
. ${KRONA_DIRE}/updateTaxonomy.sh
. ${KRONA_DIRE}/updateAccessions.sh

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

echo "Cenote-Taker2 should now run. Use \'python /path/to/Cenote-Taker2/run_cenote-taker2_200207.py\'"
#python setup.py build; python setup.py install

#Cenote-Taker2 -h