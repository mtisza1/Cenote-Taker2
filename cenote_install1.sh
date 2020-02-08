#!/bin/bash

git clone https://github.com/mtisza1/Cenote-Taker2.git

cd Cenote-Taker2
#chmod +x irf?
git clone --recursive https://github.com/deprekate/PHANOTATE.git
cd PHANOTATE; make
cd ..

conda env create --file cenote-taker2_env.yml
conda activate Cenote-Taker2
#wget hmmer DBs
#wget BLAST DBs
#wget hhsuite DBs

python setup.py build; python setup.py install

Cenote-Taker2 -h