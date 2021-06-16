#!/bin/bash

# provide exactly 1 argument. Either 'default' or an absolute path to the desired conda environment setup directory

read -p 'Hello! This installation will require 32GB of storage in your Conda environment directory and 100GB of storage in your current directory. Are you confident that you have the available space? (yes or no): ' installrvar
if [ "$installrvar" == "yes" ] ; then
	echo " "
else
	"You didn't say \"yes\". Exiting"
	exit
fi
read -p 'This installation will require 4 CPUs. Are you confident that you have 4 CUs available? (yes or no): ' cpurvar
if [ "$cpurvar" == "yes" ] ; then
	echo " "
else
	"You didn't say \"yes\". Exiting"
	exit
fi
read -p 'This installation will download and format large databases. It could take up to 2 hours. Are you ready to begin? (yes or no): ' timervar
if [ "$timervar" == "yes" ] ; then
	echo " "
	echo "Beginning install!"
	sleep 3s
else
	"You didn't say \"yes\". Exiting"
	exit
fi



if [ -s Cenote-Taker2/run_cenote-taker2.0.1.py ] ; then
	cd Cenote-Taker2
	git pull
else
	git clone https://github.com/mtisza1/Cenote-Taker2.git

	cd Cenote-Taker2
fi
#chmod +x irf307.linux.exe
# cloning PHANOTATE
if [ -s PHANOTATE/phanotate.py ] ; then
	echo "phanotate already present"
else
	git clone --recursive https://github.com/deprekate/PHANOTATE.git
	cd PHANOTATE; # git checkout version from feb 8 2020
	make
	cd ..
fi

if [ -s last-1047/README.txt ] ; then
	echo "last already present"
else
	wget http://last.cbrc.jp/last-1047.zip
	unzip last-1047.zip
	cd last-1047/
	make
	cd ..
	rm last-1047.zip
fi

eval "$(conda shell.bash hook)"
if [ "$1" == "default" ] ; then
	conda info --envs | if grep -q "cenote-taker2_env" ; then
		conda env remove --name cenote-taker2_env
	fi
	conda env create --file cenote-taker2_env.yml
	conda activate cenote-taker2_env
elif [ -d "$1" ] ; then
	conda info --envs | if grep -q "${1}\/cenote-taker2_env" ; then
		conda env remove -p ${1}/cenote-taker2_env
	fi
	conda create -p ${1}/cenote-taker2_env -c defaults -c bioconda -c AgBiome python=3.6 prodigal=2.6.3 BWA=0.7.17 samtools=1.3 mummer=3.23 circlator=1.5.5 blast=2.9.0 bioawk=1.0 entrez-direct=13.3 krona=2.7.1 hmmer=3.3 bowtie2=2.3.5 trnascan-se=2.0.5 bbtools=37.62 tbl2asn=25.7 emboss=6.6.0 cmake=3.14.0 numpy=1.18.1 pandas=1.0.0 matplotlib=3.1.3
	conda activate ${1}/cenote-taker2_env
else
	echo "no proper option given for conda environment directory. Exiting."
	exit
fi

#conda create -n cenote-taker2_env -c defaults -c bioconda -c AgBiome python=3.6 prodigal=2.6.3 BWA=0.7.17 samtools=1.3 mummer=3.23 circlator=1.5.5 blast=2.9.0 bioawk=1.0 entrez-direct=13.3 krona=2.7.1 hmmer=3.3 bowtie2=2.3.5 trnascan-se=2.0.5 bbtools=37.62 tbl2asn=25.7 emboss=6.6.0 cmake=3.14.0 numpy=1.18.1 pandas=1.0.0 matplotlib=3.1.3 

if [ "$1" == "default" ] ; then
	conda info --envs | sed 's/ \+/ /g' | if grep -q "cenote-taker2_env \*" ; then 
		echo "cenote-taker2_env loaded" ; 
	else 
		echo "cenote-taker2_env not loaded correctly" ;
		exit 
	fi
else
	conda info --envs | sed 's/ \+/ /g' | if grep -q "\* ${1}\/cenote-taker2_env" ; then 
		echo "cenote-taker2_env loaded" ; 
	else 
		echo "cenote-taker2_env not loaded correctly" ;
		exit 
	fi
fi

# installing fasthpath with pip to avoid extra channels in environment
pip install fastpath 

# getting hh-suite from github. The anaconda package of hh-suite causes conflicts
if [ -s hh-suite/README.md ] ; then
	echo "hh-suite already present"
else
	git clone https://github.com/soedinglab/hh-suite.git
	cd hh-suite # git checkout version from feb 8 2020
	mkdir -p build && cd build
	cmake -DCMAKE_INSTALL_PREFIX=. ..
	make -j 4 && make install
	export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
	cd ..
	cd ..
fi



CT2_DIR=$PWD
KRONA_DIR=$( which python | sed 's/bin\/python/opt\/krona/g' )
cd ${KRONA_DIR}
sh updateTaxonomy.sh
cd ${KRONA_DIR}
sh updateAccessions.sh
cd ${CT2_DIR}

#wget hmmer DBs
if [ -s hmmscan_DBs/useful_hmms_baits_and_not2a.h3p ] ; then
	echo "HMM databases appear to be already present"
else
	wget https://zenodo.org/record/4966268/files/hmmscan_DBs.tgz
	tar -xvf hmmscan_DBs.tgz
	rm hmmscan_DBs.tgz
fi
#wget BLAST DBs
if [ -s blast_DBs/virus_adinto_polinton_prot_190925.psq ] ; then
	echo "BLAST databases appear to be already present"
else
	wget https://zenodo.org/record/3660538/files/blast_DBs.tgz
	tar -xvf blast_DBs.tgz
	rm blast_DBs.tgz
fi
#wget hhsuite DBs
## NCBI_CD
if [ -s NCBI_CD/NCBI_CD_a3m.ffdata ] ; then
	echo "NCBI_CD hhsuite databases appear to be already present"
else	
	wget https://zenodo.org/record/3660537/files/NCBI_CD_hhsuite.tgz
	tar -xvf NCBI_CD_hhsuite.tgz
	rm NCBI_CD_hhsuite.tgz
fi
## pfam_32
if [ -s pfam_32_db/pfam_a3m.ffdata ] ; then
	echo "pfam_32 hhsuite databases appear to be already present"
else
	mkdir pfam_32_db && cd pfam_32_db
	wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_32.0.tar.gz
	tar -xvf pfamA_32.0.tar.gz
	rm pfamA_32.0.tar.gz
	cd ..
fi
## pdb70_latest
if [ -s pdb70/pdb70_a3m.ffdata ] ; then
	echo "pdb70 hhsuite databases appear to be already present"
else
	mkdir pdb70 && cd pdb70
	wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz
	tar -xvf pdb70_from_mmcif_latest.tar.gz
	rm pdb70_from_mmcif_latest.tar.gz
	cd ..
fi
#wget cdd rpsblast db
if [ -s cdd_rps_db/Cdd.rps ] ; then
	echo "RPSBLAST CDD databases appear to be already present"
else
	mkdir cdd_rps_db && cd cdd_rps_db
	wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
	tar -xvf Cdd_LE.tar.gz
	rm Cdd_LE.tar.gz
	cd ..
fi

echo "Cenote-Taker2 should now run. Use: python /path/to/Cenote-Taker2/run_cenote-taker2.py"
#python setup.py build; python setup.py install

#Cenote-Taker2 -h