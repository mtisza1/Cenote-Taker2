#!/bin/bash

# Cenote-Taker Logo
echo "$(tput setaf 2)00000000000000000000000000"
echo "00000000000000000000000000"
echo "000000000$(tput setaf 4)^^^^^^^^$(tput setaf 2)000000000"
echo "000000$(tput setaf 4)^^^^^^^^^^^^^^$(tput setaf 2)000000"
echo "00000$(tput setaf 4)^^^^^$(tput setaf 3)CENOTE$(tput setaf 4)^^^^^$(tput setaf 2)00000"
echo "00000$(tput setaf 4)^^^^^$(tput setaf 3)TAKER!$(tput setaf 4)^^^^^$(tput setaf 2)00000"
echo "00000$(tput setaf 4)^^^^^^^^^^^^^^^^$(tput setaf 2)00000"
echo "000000$(tput setaf 4)^^^^^^^^^^^^^^$(tput setaf 2)000000"
echo "000000000$(tput setaf 4)^^^^^^^^$(tput setaf 2)000000000"
echo "00000000000000000000000000"
echo "00000000000000000000000000$(tput sgr 0)"

echo " "
echo "Version 2.1.5"
echo " "
echo "Fun fact: As of this version, the$(tput setaf 4) virion $(tput sgr 0)database is used by default. Wow!"
sleep 2s


# Setting input parameters  
original_contigs=$1
F_READS=$2
R_READS=$3
run_title=$4
isolation_source=$5
collection_date=$6
metagenome_type=$7
srr_number=$8
srx_number=$9
biosample=${10}
bioproject=${11}
template_file=${12}
circ_length_cutoff=${13}
linear_length_cutoff=${14}
virus_domain_db=${15}
LIN_MINIMUM_DOMAINS=${16}
handle_knowns=${17}
ASSEMBLER=${18}
MOLECULE_TYPE=${19}
HHSUITE_TOOL=${20}
DATA_SOURCE=${21}
BLASTP=${22}
PROPHAGE=${23}
FOR_PLASMIDS=${24}
BLASTN_DB=${25}
CENOTE_SCRIPT_DIR=${26}
CIRC_MINIMUM_DOMAINS=${27}
SCRATCH_DIR=${28}
MEM=${29}
CPU=${30}
ENFORCE_START_CODON=${31}
ORF_WITHIN=${32}
ANNOTATION_MODE=${33}
CRISPR_FILE=${34}
base_directory=$PWD
CENOTE_DBS=${35}
WRAP=${36}
HALLMARK_TAX=${37}
EXACT=${38}

if [ "$ANNOTATION_MODE" == "True" ] ; then
	LIN_MINIMUM_DOMAINS=0
	CIRC_MINIMUM_DOMAINS=0
	circ_length_cutoff=1000
	linear_length_cutoff=1
	PROPHAGE="False"
fi

#># check for correct format run name

echo "@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "Your specified arguments:"
echo "original contigs:                  $original_contigs"
echo "forward reads:                     $F_READS"
echo "reverse reads:                     $R_READS"
echo "title of this run:                 $run_title"
echo "Isolate source:                    $isolation_source"
echo "collection date:                   $collection_date"
echo "metagenome_type:                   $metagenome_type"
echo "SRA run number:                    $srr_number"
echo "SRA experiment number:             $srx_number"
echo "SRA sample number:                 $biosample"
echo "Bioproject number:                 $bioproject"
echo "template file:                     $template_file"
echo "minimum circular contig length:    $circ_length_cutoff"
echo "minimum linear contig length:      $linear_length_cutoff"
echo "virus domain database:             $virus_domain_db"
echo "min. viral hallmarks for linear:   $LIN_MINIMUM_DOMAINS"
echo "min. viral hallmarks for circular: $CIRC_MINIMUM_DOMAINS"
echo "handle known seqs:                 $handle_knowns"
echo "contig assembler:                  $ASSEMBLER"
echo "DNA or RNA:                        $MOLECULE_TYPE"
echo "HHsuite tool:                      $HHSUITE_TOOL"
echo "original or TPA:                   $DATA_SOURCE"
echo "Do BLASTP?:                        $BLASTP"
echo "Do Prophage Pruning?:              $PROPHAGE"  
echo "Filter out plasmids?:              $FOR_PLASMIDS"
echo "Run BLASTN against nt?             $BLASTN_DB"
echo "Location of Cenote scripts:        $CENOTE_SCRIPT_DIR"
echo "Location of scratch directory:     $SCRATCH_DIR"
echo "GB of memory:                      $MEM"
echo "number of CPUs available for run:  $CPU"
echo "Annotation only mode?              $ANNOTATION_MODE"
echo "Cenote-Taker 2 DBs directory:      $CENOTE_DBS"
echo "Wrap circular contigs?:            $WRAP"
echo "Taxonomy for each hallmark?:       $HALLMARK_TAX"
echo "Exact match DTRs?:                 $EXACT"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@"

#checking validity of run_title
if [[ "$run_title" =~ ^[a-zA-Z0-9_]+$ ]] && [ ${#run_title} -le 18 ] ; then 
	echo $run_title ; 
else
	echo "$run_title is not a valid name for the run title ( -r argument)"
	echo " the run title needs to be only letters, numbers and underscores (_) and 18 characters or less. Exiting."
	exit
fi

if [ "${SCRATCH_DIR}" == "none" ] ; then
	echo "scratch space will not be used in this run"
	if [ -s ${CENOTE_DBS}/NCBI_CD/NCBI_CD_a3m.ffdata ] ; then
		CD_HHSUITE="${CENOTE_DBS}/NCBI_CD/NCBI_CD"
	else
		CD_HHSUITE=""
	fi
	if [ -s ${CENOTE_DBS}/pfam_32_db/pfam_a3m.ffdata ] ; then
		PFAM_HHSUITE="${CENOTE_DBS}/pfam_32_db/pfam"
	else
		PFAM_HHSUITE=""
	fi
	if [ -s ${CENOTE_DBS}/pdb70/pdb70_a3m.ffdata ] ; then
		PDB_HHSUITE="${CENOTE_DBS}/pdb70/pdb70"
	else
		PDB_HHSUITE=""
	fi

	echo "HHsuite database locations:"
	echo $CD_HHSUITE
	echo $PFAM_HHSUITE
	echo $PDB_HHSUITE

# downloading taxdump database if it doesn't exist
if [ ! -s ${CENOTE_DBS}/taxdump/names.dmp ] ; then
	echo "the required taxdump file (new requirement as of Cenote-Taker 2.1.5) wasn't found. It will be downloaded and upacked at ${CENOTE_DBS}/taxdump/"
	if [ ! -d ${CENOTE_DBS}/taxdump ] ; then
		mkdir ${CENOTE_DBS}/taxdump
	fi
	cd ${CENOTE_DBS}/taxdump
	wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	tar -xf taxdump.tar.gz
	cd $base_directory
fi

# make sure all required database directories exist
if [ ! -d ${CENOTE_DBS}/blast_DBs/ ] ; then
	echo "Can not find blast databases at ${CENOTE_DBS}/blast_DBs/"
	echo "Exiting"
	exit
fi
if [ ! -d ${CENOTE_DBS}/hmmscan_DBs/ ] ; then
	echo "Can not find HMM databases at ${CENOTE_DBS}/hmmscan_DBs/"
	echo "Exiting"
	exit
fi
if [ ! -d ${CENOTE_DBS}/cdd_rps_db/ ] ; then
	echo "Can not find RPSBLAST databases at ${CENOTE_DBS}/hmmscan_DBs/"
	echo "Exiting"
	exit
fi
if [ ! -d ${CENOTE_DBS}/taxdump/ ] ; then
	echo "Can not find Taxonomy databases at ${CENOTE_DBS}/hmmscan_DBs/"
	echo "Exiting"
	exit
fi


elif [ -d ${SCRATCH_DIR}/ ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: setting up lscratch databases: " $MDYT
	if [ ! -s ${SCRATCH_DIR}/NCBI_CD/NCBI_CD_a3m.ffdata ] ; then
		mkdir ${SCRATCH_DIR}/NCBI_CD
		cp ${CENOTE_DBS}/NCBI_CD/NCBI_CD* ${SCRATCH_DIR}/NCBI_CD/
	fi
	if [ ! -s ${SCRATCH_DIR}/pfam_32_db/pfam_a3m.ffdata ] ; then	
		mkdir ${SCRATCH_DIR}/pfam_32_db
		cp ${CENOTE_DBS}/pfam_32_db/pfam* ${SCRATCH_DIR}/pfam_32_db/
	fi
	if [ ! -s ${SCRATCH_DIR}/pdb70/pdb70_a3m.ffdata ] ; then		
		mkdir ${SCRATCH_DIR}/pdb70
		cp ${CENOTE_DBS}/pdb70/pdb70* ${SCRATCH_DIR}/pdb70/
	fi
	CD_HHSUITE="${SCRATCH_DIR}/NCBI_CD/NCBI_CD"
	PFAM_HHSUITE="${SCRATCH_DIR}/pfam_32_db/pfam"
	PDB_HHSUITE="${SCRATCH_DIR}/pdb70/pdb70"
#	mkdir /lscratch/$SLURM_JOB_ID/viral
#	cp /fdb/blastdb/viral.* /lscratch/$SLURM_JOB_ID/viral/
else
	echo "CAN'T FIND SCRATCH FOLDER. Specify correct folder or don't use scratch space. Exiting."
	exit
fi

# setting HHSUITE database argument string
HHSUITE_DB_STR=""
if [ -n "$CD_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${CD_HHSUITE} "
fi
if [ -n "$PFAM_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${PFAM_HHSUITE} "
fi
if [ -n "$PDB_HHSUITE" ] ; then
	HHSUITE_DB_STR="${HHSUITE_DB_STR}-d ${PDB_HHSUITE}"
fi
if [ ! -n "$HHSUITE_DB_STR" ] ; then
	echo "HHsuite databases not found at ${CENOTE_DBS}"
	echo "$HHSUITE_TOOL will not be run"
	HHSUITE_TOOL="none"
fi

# looking for template file and contigs in working directory, or else copying them there
if [ -s ${base_directory}/${template_file} ] ; then 
	echo ${base_directory}/${template_file} ; 
else  
	cp ${template_file} ${base_directory}/ ;
	template_file=$( basename $template_file ) 
fi

if [ -s ${base_directory}/${original_contigs} ] ; then 
	echo ${base_directory}/${original_contigs} ; 
else  
	cp ${original_contigs} ${base_directory}/ ; 
	original_contigs=$( basename $original_contigs ) 
fi


if [ -s ${base_directory}/${CRISPR_FILE} ] ; then 
	echo ${base_directory}/${CRISPR_FILE} ; 
elif [ "${CRISPR_FILE}" == "none" ] ; then
	echo "no CRISPR file given"
else  
	cp ${CRISPR_FILE} ${base_directory}/ ; 
	CRISPR_FILE=$( basename $CRISPR_FILE ) 
fi

if [ "$PROPHAGE" == "True" ] && [ $LIN_MINIMUM_DOMAINS -le 0 ] ; then
	echo "Prophage pruning requires --lin_minimum_hallmark_genes >= 1. changing to:"
	echo "--lin_minimum_hallmark_genes 1"
	LIN_MINIMUM_DOMAINS=1
fi


MDYT=$( date +"%m-%d-%y---%T" )
echo "time update: locating inputs: " $MDYT

# looking for template file and contigs in working directory, or else copying them there


if [ -s ${base_directory}/${original_contigs} ] ; then 
	echo ${base_directory}/${original_contigs} ; 
else  
	cp ${original_contigs} ${base_directory}/ ; 
	original_contigs=$( basename $original_contigs ) 

fi

# Making output folder
if [ ! -d "$run_title" ]; then
	mkdir "$run_title"
else
	rand_dir=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )
	DAY1=$( date +"%m-%d-%y" )
	mv ${run_title}/ ${run_title}_old_${DAY1}_${rand_dir} 
	mkdir "$run_title"
fi 

if [ ${original_contigs: -3} == ".fa" ]; then
	echo "renaming $original_contigs to ${original_contigs}sta"
	mv $original_contigs ${original_contigs}sta
	original_contigs=${original_contigs}sta
fi
if [ ${original_contigs: -4} == ".fna" ]; then
	echo "renaming $original_contigs to ${original_contigs%fna}fasta"
	mv $original_contigs ${original_contigs%fna}fasta
	original_contigs=${original_contigs%fna}fasta
fi
if [ ${original_contigs: -4} == ".fsa" ]; then
	echo "renaming $original_contigs to ${original_contigs%fsa}fasta"
	mv $original_contigs ${original_contigs%fsa}fasta
	original_contigs=${original_contigs%fsa}fasta
fi
# Removing contigs under $circ_length_cutoff nts and detecting circular contigs
if [ $circ_length_cutoff -gt $linear_length_cutoff ] ; then
	LENGTH_MINIMUM=$linear_length_cutoff
else
	LENGTH_MINIMUM=$circ_length_cutoff

fi

#LASTDBQ=$( find $CENOTE_SCRIPT_DIR -type f -wholename "$CENOTE_SCRIPT_DIR/last-*/lastdb" )
#LASTALQ=$( find $CENOTE_SCRIPT_DIR -type f -wholename "$CENOTE_SCRIPT_DIR/last-*/lastal" )
if [ ${original_contigs: -6} == ".fasta" ]; then
	echo "$(tput setaf 5)File with .fasta extension detected, attempting to keep contigs over $LENGTH_MINIMUM nt and find circular sequences with apc.pl$(tput sgr 0)"
	bioawk -v run_var="$run_title" -v contig_cutoff="$LENGTH_MINIMUM" -c fastx '{ if(length($seq) > contig_cutoff) { print ">"run_var NR" "$name; print $seq }}' $original_contigs > ${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta ;
	cd $run_title
	echo "cenote_shortcut" > ${run_title}_CONTIG_SUMMARY.tsv
	if [ "$EXACT" == "True" ] ; then
		perl ${CENOTE_SCRIPT_DIR}/apc_exact1.pl -b $run_title -c lastdb -d lastal ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta >/dev/null 2>&1
	else
		perl ${CENOTE_SCRIPT_DIR}/apc_cenote1.pl -b $run_title -c lastdb -d lastal ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta >/dev/null 2>&1
	fi
	find . -type f -name "apc_aln*" -exec rm -f {} \;
	APC_CIRCS=$( find . -maxdepth 1 -type f -name "${run_title}*.fa" )
	if [ -n "$APC_CIRCS" ] ;then
		for fa1 in $APC_CIRCS ; do
			#-#-# adding wrap option to clip DTRs or not
			if [ "$WRAP" == "True" ] ; then
				CIRC_NEW_NAME=$( head -n1 $fa1 | sed 's/|.*//g ; s/>//g ; s/ .*//g' ) ; 
				sed 's/|.*//g ; /^$/d' $fa1 > ${CIRC_NEW_NAME}.fasta
			else
				CIRC_SEQ_NAME=$( head -n1 $fa1 | sed 's/|.*//g' ) ; 
				CIRC_NEW_NAME=$( echo "$CIRC_SEQ_NAME" | sed 's/>//g ; s/ .*//g' )
				grep -A1 "^$CIRC_SEQ_NAME" ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | sed '/--/d' > ${CIRC_NEW_NAME}.fasta
			fi
			echo "${CIRC_NEW_NAME}.fasta has DTRs/circularity"
			rm -f $fa1
		done 
	else
		echo "No circular contigs detected."
	fi		
elif [ ${original_contigs: -6} == ".fastg" ]; then
	bioawk -v contig_cutoff="$LENGTH_MINIMUM" -c fastx '{ if(length($seq) > contig_cutoff) {print }}' $original_contigs | grep "[a-zA-Z0-9]:\|[a-zA-Z0-9];" | grep -v "':" | awk '{ print ">"$1 ; print $2 }' | sed 's/:.*//g; s/;.*//g' | bioawk -v run_var="$run_title" -c fastx '{ print ">"run_var NR" "$name; print $seq }' > ${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta
	cd $run_title
	echo "cenote_shortcut" > ${run_title}_CONTIG_SUMMARY.tsv
	perl ${CENOTE_SCRIPT_DIR}/apc_cenote1.pl -b $run_title -c $LASTDBQ -d $LASTALQ ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta >/dev/null 2>&1
	rm -f apc_aln*
	APC_CIRCS=$( find . -maxdepth 1 * -type f -name "${run_title}*.fa" )
	if [ -n "$APC_CIRCS" ] ;then
		for fa1 in $APC_CIRCS ; do
			if [ "$WRAP" == "True" ] ; then
				CIRC_NEW_NAME=$( head -n1 $fa1 | sed 's/|.*//g ; s/>//g ; s/ .*//g' ) ; 
				sed 's/|.*//g ; /^$/d' $fa1 > ${CIRC_NEW_NAME}.fasta
			else
				CIRC_SEQ_NAME=$( head -n1 $fa1 | sed 's/|.*//g' ) ; 
				CIRC_NEW_NAME=$( echo "$CIRC_SEQ_NAME" | sed 's/>//g ; s/ .*//g' )
				grep -A1 "^$CIRC_SEQ_NAME" ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | sed '/--/d' > ${CIRC_NEW_NAME}.fasta
			fi
			echo "${CIRC_NEW_NAME}.fasta has DTRs/circularity"
			rm -f $fa1
		done 
	else
		echo "No circular contigs detected."
	fi	
else
	echo "$(tput setaf 4)File with .fasta of .fastg extension not detected as first input. Exiting.$(tput sgr 0)" ;
	exit
fi

# Removing cirles that are smaller than user specified cutoff
CIRC_CONTIGS=$( find . -maxdepth 1 -type f -name "*.fasta" )
if [ ! -z "$CIRC_CONTIGS" ] ;then
	for CIRCLE1 in $CIRC_CONTIGS ; do
		CIRCLE1_LENGTH=$( bioawk -c fastx '{print length($seq) }' $CIRCLE1 )
		if [[ $CIRCLE1_LENGTH -lt $circ_length_cutoff ]] ; then
			mv $CIRCLE1 ${CIRCLE1%.fasta}.too_short.fasta 
		fi
	done
fi
rm -f *.too_short.fasta

# Aligning reads to contigs
FIRST_F_READS=$( echo "$F_READS" | sed 's/,.*//g' )
if [ ! -s "$FIRST_F_READS" ] ; then
	echo "no reads provided or reads not found"
else
	echo "$(tput setaf 4)Aligning provided reads to contigs over cutoff to determine coverage. $(tput sgr 0)" 
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: making bowtie2 indices " $MDYT
	mkdir bt2_indices ; 
	#### update the script to allow for contigs from other directories
	bowtie2-build ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta bt2_indices/${run_title}_bt2_index >/dev/null 2>&1
	echo "$(tput setaf 4)Aligning reads to BowTie2 index. $(tput sgr 0)" 
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: aligning reads to bowtie2 indices " $MDYT
	if [ ! -s "$R_READS" ] ;then
		bowtie2 -q -p${CPU} -x bt2_indices/${run_title}_bt2_index -U $F_READS -S reads_to_all_contigs_over${LENGTH_MINIMUM}nt.sam --very-fast
	else
		bowtie2 -q -p${CPU} -x bt2_indices/${run_title}_bt2_index -1 $F_READS -2 $R_READS -S reads_to_all_contigs_over${LENGTH_MINIMUM}nt.sam --very-fast
	fi
	echo "$(tput setaf 4)Calculating read coverage of each contig with BBTools Pileup. $(tput sgr 0)" 
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running BBTools Pileup " $MDYT
	pileup.sh in=reads_to_all_contigs_over${LENGTH_MINIMUM}nt.sam out=reads_to_all_contigs_over${LENGTH_MINIMUM}nt.coverage.txt
	rm -f reads_to_all_contigs_over${LENGTH_MINIMUM}nt.sam
fi

# Detecting whether any circular contigs were present
original_fastas=$( find . -maxdepth 1 -type f -name "*.fasta" )
# "$(tput setaf 5)$var1$(tput sgr 0)"

if [ -z "$original_fastas" ] ; then
	echo "$(tput setaf 4)No circular fasta files detected. $(tput sgr 0)" 
	#exit
	mkdir other_contigs	
	if [ ${original_contigs: -6} == ".fasta" ]; then
		grep -A1 "^>" ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | sed 's/--//g' | bioawk -v contig_cutoff="$linear_length_cutoff" -c fastx '{ if (length($seq) > contig_cutoff) { print ">"$name" "$4 ; print $seq }}' > other_contigs/all_non_circular.fasta
	elif [ ${original_contigs: -6} == ".fastg" ]; then
		grep -A1 "^>" ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta | sed 's/--//g' | bioawk -v contig_cutoff="$linear_length_cutoff" -c fastx '{ if (length($seq) > contig_cutoff) { print ">"$name" "$4 ; print $seq }}' > other_contigs/all_non_circular.fasta
	fi			
	if [ -s other_contigs/all_non_circular.fasta ] ; then
		grep "^>" other_contigs/all_non_circular.fasta | sed 's/>//g' | cut -d " " -f1 | while read LINE ; do 
			grep -A1 "$LINE [^\s]" other_contigs/all_non_circular.fasta > other_contigs/$LINE.fasta ; 
		done
	fi
else
	echo "$(tput setaf 5)Circular fasta file(s) detected$(tput sgr 0)"
	echo " "
	cat *fasta > all_circular_contigs_${run_title}.fna
# Putting non-circular contigs in a separate directory
	echo "$(tput setaf 4)Putting non-circular contigs in a separate directory $(tput sgr 0)" 

	mkdir other_contigs

	for CIRCLE in $original_fastas ; do
		grep "^>" $CIRCLE | sed 's/|.*//g' >> circular_contigs_spades_names.txt
	done
	if [ ${original_contigs: -6} == ".fasta" ]; then
		LINEAR_COUNT=$( grep -v -f circular_contigs_spades_names.txt ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | grep -A1 "^>" | wc -l )
		if [ $LINEAR_COUNT -ge 1 ] ; then
			grep -v -f circular_contigs_spades_names.txt ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | grep -A1 "^>" | sed 's/--//g' | bioawk -v contig_cutoff="$linear_length_cutoff" -c fastx '{ if (length($seq) > contig_cutoff) { print ">"$name" "$4 ; print $seq }}' > other_contigs/all_non_circular.fasta
		else
			echo "No linear contigs over ${LENGTH_MINIMUM}nt found in run"
		fi
	elif [ ${original_contigs: -6} == ".fastg" ]; then
		LINEAR_COUNT=$( grep -v -f circular_contigs_spades_names.txt ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta | grep -A1 "^>" | wc -l )
		if [ $LINEAR_COUNT -ge 1 ] ; then		
			grep -v -f circular_contigs_spades_names.txt ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta | grep -A1 "^>" | sed 's/--//g' | bioawk -v contig_cutoff="$linear_length_cutoff" -c fastx '{ if (length($seq) > contig_cutoff) { print ">"$name" "$4 ; print $seq }}' > other_contigs/all_non_circular.fasta
		else
			echo "No linear contigs over ${LENGTH_MINIMUM}nt found in run"
		fi
	fi
	if [ -s other_contigs/all_non_circular.fasta ] ; then
		grep "^>" other_contigs/all_non_circular.fasta | sed 's/>//g' | cut -d " " -f1 | while read LINE ; do 
			grep -A1 "$LINE [^\s]" other_contigs/all_non_circular.fasta > other_contigs/$LINE.fasta ; 
		done
	fi

fi

cd other_contigs
CONTIGS_NON_CIRCULAR=$( find . -maxdepth 1 -type f -name "*[0-9].fasta" )

if [ ! -z "$CONTIGS_NON_CIRCULAR" ] ;then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running IRF for ITRs in non-circular contigs" $MDYT
	for NONCIR in $CONTIGS_NON_CIRCULAR ; do
		LEN_CHECKQ=$( cat $NONCIR | bioawk -c fastx '{ if(length($seq) > 4000) { print $name }}' ) ; 
		if [ ! -z "$LEN_CHECKQ" ] ; then
			${CENOTE_SCRIPT_DIR}/irf307.linux.exe $NONCIR 2 3 5 80 10 40 500000 10000 -d -h >/dev/null 2>&1
		fi
	done
	mkdir ../ITR_containing_contigs
	DAT_FILES=$( find  * -maxdepth 0 -type f -name "*.dat" )
	if [ -n "$DAT_FILES" ] ; then
		for DAT in $DAT_FILES ; do 

			if grep -q "^[0-9]" $DAT ; then 
				LENGTH=$( bioawk -c fastx '{ print length($seq) }' ${DAT%.2.3.5.80.10.40.500000.10000.dat} | sed 's/ //g; s/\///g' ) ;
				ITR_STARTS=$( grep "^[0-9]" $DAT | cut -d " " -f1 ) ; 
				ITR_ENDS=$( grep "^[0-9]" $DAT | cut -d " " -f5 ) ; 
				LOW_START=$( echo $ITR_STARTS | tr " " "\n" | sort -g | head -n1| sed 's/ //g' ) ; 
				HIGH_END=$( echo $ITR_ENDS | tr " " "\n" | sort -g | tail -n1 | sed 's/ //g' ) ; 
				if [ $LOW_START -gt 0 ] && [ $HIGH_END -gt 0 ] ; then 
					HIGH_DIST=$(( $LENGTH-$HIGH_END )) ; 
					if [ $LOW_START -lt 1000 ] && [ $HIGH_DIST -lt 1000 ] ; then 
						echo ${DAT%.2.3.5.80.10.40.500000.10000.dat} "contains ITRs" ; 
						#echo $LENGTH "5-prime ITR:" $LOW_START "3-prime ITR:" $HIGH_END ; 
						mv ${DAT%.2.3.5.80.10.40.500000.10000.dat} ../ITR_containing_contigs/${DAT%.fasta.2.3.5.80.10.40.500000.10000.dat}.fasta
						#echo "$(tput setaf 4) Making ITR .tbl file $(tput sgr 0)"
						L_END_A=$( grep "^$LOW_START" $DAT | cut -d " " -f2 )
						L_START_B=$( grep "^$LOW_START" $DAT | cut -d " " -f4 )
						L_END_B=$( grep "^$LOW_START" $DAT | cut -d " " -f5 )
						H_START_A=$( grep " $HIGH_END " $DAT | cut -d " " -f1 )
						H_END_A=$( grep " $HIGH_END " $DAT | cut -d " " -f2 )		
						H_START_B=$( grep " $HIGH_END " $DAT | cut -d " " -f4 )
						echo -e "$LOW_START\t""$L_END_A\t""repeat_region\n""\t\t\trpt_type\tITR\n""$L_START_B\t""$L_END_B\t""repeat_region\n""\t\t\trpt_type\tITR\n""$H_START_A\t""$H_END_A\t""repeat_region\n""\t\t\trpt_type\tITR\n""$H_START_B\t""$HIGH_END\t""repeat_region\n""\t\t\trpt_type\tITR\n" >> ../ITR_containing_contigs/${DAT%.fasta.2.3.5.80.10.40.500000.10000.dat}.ITR.tbl; 
					fi ; 
				fi ; 
			fi ; 
		done
	fi
	CONTIGS_NON_CIRCULAR=$( find . -maxdepth 1 -type f -name "*[0-9].fasta" )
	mkdir ../no_end_contigs_with_viral_domain
	if [ $LIN_MINIMUM_DOMAINS -le 0 ] ; then
		for LIN in $CONTIGS_NON_CIRCULAR ; do
			mv ${LIN} ../no_end_contigs_with_viral_domain/${LIN%.fasta}.fna
		done
	else
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running prodigal on linear contigs " $MDYT
		echo "$CONTIGS_NON_CIRCULAR" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t prodigal -a {}.AA.fasta -i {}.fasta -p meta -q >/dev/null 2>&1
		for NO_END in $CONTIGS_NON_CIRCULAR ; do 
			sed 's/ /@/g' ${NO_END%.fasta}.AA.fasta | bioawk -c fastx '{print}' | while read LINE ; do 
				START_BASE=$( echo "$LINE" | cut -d "#" -f 2 | sed 's/@//g' ) ; 
				END_BASE=$( echo "$LINE" | cut -d "#" -f 3 | sed 's/@//g' ) ; 
				ORF_NAME=$( echo "$LINE" | cut -d "#" -f 1 | sed 's/@//g; s/\./_/g' ) ; 
				AA_SEQ=$( echo "$LINE" | cut -f2 ) ;
				ORIG_CONTIG=$( grep ">" ${NO_END} | cut -d " " -f2 )
				if echo $AA_SEQ | grep -q "\*" ; then
					INC3=""
				else
					INC3="3primeInc"
				fi
				FAA=${AA_SEQ:0:1}
				if [ "$FAA" != "M" ] && [ $START_BASE -le 3 ] ; then
					INC5="5primeInc"
				else
					INC5=""
				fi
				echo ">"${ORF_NAME} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${NO_END%.fasta}.AA.sorted.fasta
		done
		LIN_NO_ITR_AA=$( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" )
		if [ -n "$LIN_NO_ITR_AA" ] ; then
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running linear contigs with hmmscan against virus hallmark gene database: $virus_domain_db " $MDYT	
			#cat $( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" ) > all_large_genome_proteins.AA.fasta
			for LIN in $LIN_NO_ITR_AA ; do
				cat $LIN
			done > all_large_genome_proteins.AA.fasta
			TOTAL_AA_SEQS=$( grep -F ">" all_large_genome_proteins.AA.fasta | wc -l | bc )
			if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
				AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
				if [ $AA_SEQS_PER_FILE = 0 ] ; then
					AA_SEQS_PER_FILE=1
				fi
				awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_LARGE_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_large_genome_proteins.AA.fasta
				SPLIT_AA_LARGE=$( find . -maxdepth 1 -type f -name "SPLIT_LARGE_GENOME_AA_*.fasta" )
				if  [[ $virus_domain_db = "standard" ]] ; then
					echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
					echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
				elif [[ $virus_domain_db = "virion" ]]; then
					echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
				elif [[ $virus_domain_db = "rna_virus" ]]; then
					echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
				else
					echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
					rm -f ./*{0..9}.fasta
					break
				fi		
				HMM_REP_NUMEBR=$( find . -maxdepth 1 -type f -name "SPLIT_LARGE_GENOME_AA_*AA.hmmscan_replicate.out" | wc -l )
				if [[ $FOR_PLASMIDS = "True" ]]; then
					if [ $HMM_REP_NUMEBR -gt 0 ] ; then
						cat SPLIT_LARGE_GENOME_AA_*AA.hmmscan.out SPLIT_LARGE_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > LARGE_GENOME_COMBINED.AA.hmmscan.sort.out
					else
						cat SPLIT_LARGE_GENOME_AA_*AA.hmmscan.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > LARGE_GENOME_COMBINED.AA.hmmscan.sort.out
					fi
				else
					if [ $HMM_REP_NUMEBR -gt 0 ] ; then
						cat SPLIT_LARGE_GENOME_AA_*AA.hmmscan.out SPLIT_LARGE_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > LARGE_GENOME_COMBINED.AA.hmmscan.sort.out
					else
						cat SPLIT_LARGE_GENOME_AA_*AA.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > LARGE_GENOME_COMBINED.AA.hmmscan.sort.out
					fi
				fi		
				if [ -e LARGE_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
					cut -f3 LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
						HALL_COUNT=$( grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | wc -l | bc )
						if [ $HALL_COUNT -ge $LIN_MINIMUM_DOMAINS ] ; then 
							mv "${HIT}.fasta" "../no_end_contigs_with_viral_domain/${HIT}.fna"
							grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out > ../no_end_contigs_with_viral_domain/${HIT}.AA.hmmscan.sort.out
							# ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.no_hmmscan1.fasta
							grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan.txt
							grep -v -f ${HIT}.AA.called_hmmscan.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ../no_end_contigs_with_viral_domain/${HIT}.AA.no_hmmscan1.fasta
							mv "${HIT}.AA.sorted.fasta" ../no_end_contigs_with_viral_domain/

						else
							cat ${HIT}.fasta >> non_viral_domains_contigs.fna
							rm -f ${HIT}.fasta
						fi
					done
				fi
			else
				echo "no AA seqs found in linear contigs, discover virus module"
			fi
			for REMAINDER in $CONTIGS_NON_CIRCULAR ; do
				if [ -s $REMAINDER ] ; then
					if [ "$ANNOTATION_MODE" == "True" ] ; then
						mv $REMAINDER ../no_end_contigs_with_viral_domain/${REMAINDER%.fasta}.fna
					else
						cat $REMAINDER >> non_viral_domains_contigs.fna
						rm -f $REMAINDER
					fi
				fi
			done
		else
			echo "no linear seqs without ITRs"
		fi
	fi
fi


cd ${base_directory}/${run_title}


DTR_SEQS=$( find . -maxdepth 1 -type f -regextype sed -regex "./${run_title}[0-9]\{1,6\}.fasta" | sed 's/\.\///g' )


if [ ! -z "$DTR_SEQS" ] ; then
	mkdir DTR_contigs_with_viral_domain
	if [ $CIRC_MINIMUM_DOMAINS -le 0 ] ; then

		for CIRC in $DTR_SEQS ; do
			mv ${CIRC} DTR_contigs_with_viral_domain/${CIRC%.fasta}.fna
			mv ${CIRC%.fasta}.DTR.tbl DTR_contigs_with_viral_domain/
		done
	else
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: Calling ORFs for circular/DTR sequences with prodigal " $MDYT
		echo "$DTR_SEQS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t prodigal -a {}.AA.fasta -i {}.fasta -p meta -q >/dev/null 2>&1
		for CIRC in $DTR_SEQS ; do 
			mv ${CIRC%.fasta}.DTR.tbl DTR_contigs_with_viral_domain/
			sed 's/ /@/g' ${CIRC%.fasta}.AA.fasta | bioawk -c fastx '{print}' | while read LINE ; do 
				START_BASE=$( echo "$LINE" | cut -d "#" -f 2 | sed 's/@//g' ) ; 
				END_BASE=$( echo "$LINE" | cut -d "#" -f 3 | sed 's/@//g' ) ; 
				ORF_NAME=$( echo "$LINE" | cut -d "#" -f 1 | sed 's/@//g; s/\./_/g' ) ; 
				AA_SEQ=$( echo "$LINE" | cut -f2 ) ;
				ORIG_CONTIG=$( grep ">" ${CIRC} | cut -d " " -f2 )
				if echo $AA_SEQ | grep -q "\*" ; then
					INC3=""
				else
					INC3="3primeInc"
				fi
				FAA=${AA_SEQ:0:1}
				if [ "$FAA" != "M" ] && [ $START_BASE -le 3 ] ; then
					INC5="5primeInc"
				else
					INC5=""
				fi
				echo ">"${ORF_NAME} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${CIRC%.fasta}.AA.sorted.fasta
		done	
		#cat $( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" ) > all_circular_genome_proteins.AA.fasta
		DTR_SEQS_SORT_AA=$( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" )
		for DTRQ in $DTR_SEQS_SORT_AA ; do
			cat $DTRQ
		done > all_circular_genome_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_circular_genome_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_CIRCULAR_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_circular_genome_proteins.AA.fasta
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running hmmscan on circular/DTR contigs " $MDYT
			SPLIT_AA_CIRC=$( find . -maxdepth 1 -type f -name "SPLIT_CIRCULAR_AA_*.fasta" )
			if  [[ $virus_domain_db = "standard" ]] ; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_replication_clusters3 {}.fasta	>/dev/null 2>&1
			elif [[ $virus_domain_db = "virion" ]]; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta	>/dev/null 2>&1	
			elif [[ $virus_domain_db = "rna_virus" ]]; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
			else
				echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
				rm -f ./*{0..9}.fasta
				break
			fi
			HMM_REP_NUMEBR=$( find . -maxdepth 1 -type f -name "SPLIT_CIRCULAR_AA_*AA.hmmscan_replicate.out" | wc -l )
			if [[ $FOR_PLASMIDS = "True" ]]; then
				if [ $HMM_REP_NUMEBR -gt 0 ] ; then
					cat SPLIT_CIRCULAR_AA_*AA.hmmscan.out SPLIT_CIRCULAR_AA_*AA.hmmscan_replicate.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out
				else
					cat SPLIT_CIRCULAR_AA_*AA.hmmscan.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out
				fi
			else
				if [ $HMM_REP_NUMEBR -gt 0 ] ; then
					cat SPLIT_CIRCULAR_AA_*AA.hmmscan.out SPLIT_CIRCULAR_AA_*AA.hmmscan_replicate.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out
				else
					cat SPLIT_CIRCULAR_AA_*AA.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out
				fi
			fi		
			if [ -s CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
				cut -f3 CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
					HALL_COUNT=$( grep "${HIT}_" CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out | wc -l | bc )
					if [ $HALL_COUNT -ge $CIRC_MINIMUM_DOMAINS ] ; then 
						mv ${HIT}.fasta DTR_contigs_with_viral_domain/${HIT}.fna
						grep "${HIT}_" CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out > DTR_contigs_with_viral_domain/${HIT}.AA.hmmscan.sort.out
						grep "${HIT}_" CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /'> ${HIT}.AA.called_hmmscan.txt
						grep -v -f ${HIT}.AA.called_hmmscan.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > DTR_contigs_with_viral_domain/${HIT}.AA.no_hmmscan1.fasta
						mv ${HIT}.AA.sorted.fasta DTR_contigs_with_viral_domain/
						mv ${HIT}.AA.fasta DTR_contigs_with_viral_domain/

					elif [ -s ${HIT}.fasta ] ; then
						sed 's/ /#/g' ${HIT}.fasta | bioawk -c fastx '{print ">"$name"#DTRs" ; print $seq}' | sed 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
						rm -f ${HIT}.fasta
						rm -f ${HIT}.AA.fasta
						rm -f ${HIT}.AA.sorted.fasta
					fi
				done
			fi
		else
			echo "No AA seqs found in circular contigs, discover viruses module"
		fi
		for REMAINDER in $DTR_SEQS ; do
			if [ -s $REMAINDER ] ; then
				if [ "$ANNOTATION_MODE" == "True" ] ; then
					mv $REMAINDER DTR_contigs_with_viral_domain/${REMAINDER%.fasta}.fna
				else
					sed 's/ /#/g' $REMAINDER | bioawk -c fastx '{print ">"$name"#DTRs" ; print $seq}' | sed 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
					rm -f ${REMAINDER%fasta}*fasta
				fi
			fi
		done
	fi
fi

ITR_SEQS=$( find * -maxdepth 1 -type f -wholename "ITR_containing_contigs/*fasta" )


if [ ! -z "$ITR_SEQS" ] ; then 
	cd ITR_containing_contigs
	ITR_SEQS=$( find . -maxdepth 1 -type f -wholename "*fasta" )
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: Calling ORFs for ITR sequences with prodigal " $MDYT
	echo "$ITR_SEQS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t prodigal -a {}.AA.fasta -i {}.fasta -p meta -q >/dev/null 2>&1
	for CIRC in $ITR_SEQS ; do 
		sed 's/ /@/g' ${CIRC%.fasta}.AA.fasta | bioawk -c fastx '{print}' | while read LINE ; do 
			START_BASE=$( echo "$LINE" | cut -d "#" -f 2 | sed 's/@//g' ) ; 
			END_BASE=$( echo "$LINE" | cut -d "#" -f 3 | sed 's/@//g' ) ; 
			ORF_NAME=$( echo "$LINE" | cut -d "#" -f 1 | sed 's/@//g; s/\./_/g' ) ; 
			AA_SEQ=$( echo "$LINE" | cut -f2 ) ;
			ORIG_CONTIG=$( grep ">" ${CIRC} | cut -d " " -f2 )
			if echo $AA_SEQ | grep -q "\*" ; then
				INC3=""
			else
				INC3="3primeInc"
			fi
			FAA=${AA_SEQ:0:1}
			if [ "$FAA" != "M" ] && [ $START_BASE -le 3 ] ; then
				INC5="5primeInc"
			else
				INC5=""
			fi
			echo ">"${ORF_NAME} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
		done > ${CIRC%.fasta}.AA.sorted.fasta
	done
	#cat $( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" ) > all_ITR_genome_proteins.AA.fasta
	ITR_SEQ_SORT_AA=$( find . -maxdepth 1 -type f -name "*.AA.sorted.fasta" )
	for ITRQ in $ITR_SEQ_SORT_AA ; do
		cat $ITRQ
	done > all_ITR_genome_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_ITR_genome_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_ITR_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_ITR_genome_proteins.AA.fasta
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running hmmscan on ITR contigs " $MDYT
		SPLIT_AA_CIRC=$( find . -maxdepth 1 -type f -name "SPLIT_ITR_AA_*.fasta" )
		if  [[ $virus_domain_db = "standard" ]] ; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
		elif [[ $virus_domain_db = "virion" ]]; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
		elif [[ $virus_domain_db = "rna_virus" ]]; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
		else
			echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
			rm -f ./*{0..9}.fasta
			break
		fi
		HMM_REP_NUMEBR=$( find . -maxdepth 1 -type f -name "SPLIT_ITR_AA_*AA.hmmscan_replicate.out" | wc -l )
		if [[ $FOR_PLASMIDS = "True" ]]; then
			if [ $HMM_REP_NUMEBR -gt 0 ] ; then
				cat SPLIT_ITR_AA_*AA.hmmscan.out SPLIT_ITR_AA_*AA.hmmscan_replicate.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > ITR_GENOME_COMBINED.AA.hmmscan.sort.out
			else
				cat SPLIT_ITR_AA_*AA.hmmscan.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > ITR_GENOME_COMBINED.AA.hmmscan.sort.out
			fi
		else
			if [ $HMM_REP_NUMEBR -gt 0 ] ; then
				cat SPLIT_ITR_AA_*AA.hmmscan.out SPLIT_ITR_AA_*AA.hmmscan_replicate.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > ITR_GENOME_COMBINED.AA.hmmscan.sort.out
			else
				cat SPLIT_ITR_AA_*AA.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > ITR_GENOME_COMBINED.AA.hmmscan.sort.out
			fi
		fi		
		if [ -e ITR_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
			cut -f3 ITR_GENOME_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
				HALL_COUNT=$( grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out | wc -l | bc )
				if [ $HALL_COUNT -ge $CIRC_MINIMUM_DOMAINS ] ; then 
					mv ${HIT}.fasta ${HIT}.fna
					grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out > ${HIT}.AA.hmmscan.sort.out
					grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan.txt
					grep -v -f ${HIT}.AA.called_hmmscan.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.AA.no_hmmscan1.fasta
					rm ${HIT}.AA.fasta

				elif [ -s ${HIT}.fasta ] ; then
					sed 's/ /#/g' ${HIT}.fasta | bioawk -c fastx '{print ">"$name"#ITRs" ; print $seq}' | sed 's/#/ /g' >> ../other_contigs/non_viral_domains_contigs.fna
					rm -f ${HIT}.fasta
				fi
			done
		fi
	else
		echo "no AA seqs found in ITR contigs, discover viruses module"
	fi
	for REMAINDER in $ITR_SEQS ; do
		if [ -s $REMAINDER ] ; then
			if [ "$ANNOTATION_MODE" == "True" ] ; then
				mv $REMAINDER ${REMAINDER%.fasta}.fna
			else
				sed 's/ /#/g' $REMAINDER | bioawk -c fastx '{print ">"$name"#ITRs" ; print $seq}' | sed 's/#/ /g' >> ../other_contigs/non_viral_domains_contigs.fna
				rm -f $REMAINDER
			fi
		fi
	done
fi


cd ${base_directory}/${run_title}



if [ -s all_circular_genome_proteins.AA.fasta ] ; then
	mv all_circular_genome_proteins.AA.fasta DTR_contigs_with_viral_domain/
fi
if [ -s CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
	mv CIRCULAR_GENOME_COMBINED.AA.hmmscan.sort.out DTR_contigs_with_viral_domain/
fi


# script for pruning no_end contigs with viral domains
LIST_OF_VIRAL_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "no_end_contigs_with_viral_domain/*fna" )
if [ ! -z "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] && [ "$PROPHAGE" == "True" ] ;then
	echo "$(tput setaf 3) Starting pruning of non-DTR/circular contigs with viral domains $(tput sgr 0)"

	. ${CENOTE_SCRIPT_DIR}/prune_linear_contigs_0.1.sh
fi

cd ${base_directory}/${run_title}

if [ -d DTR_contigs_with_viral_domain ] ; then
	cd DTR_contigs_with_viral_domain
fi

#-# rotate DTRs
rm -f *no_hmmscan1.fasta
CIRCULAR_HALLMARK_CONTIGS=$( find . -maxdepth 1 -type f -name "*fna" )
if [ -n "$CIRCULAR_HALLMARK_CONTIGS" ] ; then
	if [ "$WRAP" == "True" ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Annotating DTR contigs " $MDYT
		echo "rotating DTR contigs"
		rm *DTR.tbl
		for nucl_fa in $CIRCULAR_HALLMARK_CONTIGS ; do
			prodigal -a ${nucl_fa%.fna}.AA.fasta -i $nucl_fa -c -p meta -q >/dev/null 2>&1
			FWD_GENES=$( grep "^>" ${nucl_fa%.fna}.AA.fasta | sed 's/ # /	/g' | awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
			REV_GENES=$( grep "^>" ${nucl_fa%.fna}.AA.fasta | sed 's/ # /	/g' | awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == -1) {print $4}}' | wc -l )
			if [ $FWD_GENES -ge $REV_GENES ] && [ $FWD_GENES -ge 1 ]; then
				START_BASE=$( grep "^>" ${nucl_fa%.fna}.AA.fasta | sed 's/ # /	/g' | awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' | sort -rg -k2,2 | head -n1 | cut -f1 )
				cat $nucl_fa | seqkit restart -i ${START_BASE} > ${nucl_fa%.fna}.rotate.fasta
			elif [ $REV_GENES -ge 1 ]; then
				seqkit seq $nucl_fa --quiet -t DNA -r -p > ${nucl_fa%.fna}.rc.fna
				prodigal -a ${nucl_fa%.fna}.AA.rc.fasta -i ${nucl_fa%.fna}.rc.fna -p meta -q >/dev/null 2>&1
				RC_FWD_GENES=$( grep "^>" ${nucl_fa%.fna}.AA.rc.fasta | sed 's/ # /	/g' | awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
				if [ $RC_FWD_GENES -ge 1 ] ; then 
					START_BASE=$( grep "^>" ${nucl_fa%.fna}.AA.rc.fasta | sed 's/ # /	/g' | awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' | sort -rg -k2,2 | head -n1 | cut -f1 )
					cat ${nucl_fa%.fna}.rc.fna | seqkit restart -i ${START_BASE} > ${nucl_fa%.fna}.rotate.fasta
				else
					echo "Can't find suitable ORF to set rotation of $nucl_fa and will remain unrotated"
					cp $nucl_fa ${nucl_fa%.fna}.rotate.fasta
				fi
			else
				echo "Can't find suitable ORF to set rotation of $nucl_fa and will remain unrotated"
				cp $nucl_fa ${nucl_fa%.fna}.rotate.fasta
			fi
		done
	else
		echo "Annotating DTR contigs"
		echo "contigs will not be wrapped: --wrap False"
		echo "(rotate.fasta) files will be created for processing purposes (but these are not wrapped)"
		for nucl_fa in $CIRCULAR_HALLMARK_CONTIGS ; do
			cp $nucl_fa ${nucl_fa%.fna}.rotate.fasta
		done
	fi
fi

#-# blastx for translation decision
ROTATED_DTR_CONTIGS=$( find . -maxdepth 1 -type f -regextype sed -regex "./${run_title}[0-9]\{1,6\}.rotate.fasta" | sed 's/\.\///g' )

if [ -n "$ROTATED_DTR_CONTIGS" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running BLASTX, DTR contigs " $MDYT
	echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU -t blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -threshold 21 -word_size 5 -num_threads 1 -num_alignments 1 -db ${CENOTE_DBS}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query {}.rotate.fasta -out {}.tax_guide.blastx.out >/dev/null 2>&1
	echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta/.fasta/g' | while read nucl_fa ; do
		#-#-# Reworking taxonomy
		if [ ! -s "${nucl_fa%.fasta}.tax_guide.blastx.out" ]; then
			echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out ;
		elif grep -q "virophage" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Virophage	Unclassified Taxon" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		elif grep -q "adinto" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Adintovirus	Unclassified Taxon" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		elif grep -i -q "polinton" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Polinton-like virus	Unclassified Taxon" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		else
			ORGANISM_H=$( head -n1 ${nucl_fa%.fasta}.tax_guide.blastx.out | sed 's/\[/&\n/;s/.*\n//;s/\]/\n&/;s/\n.*//' )
			if grep -q "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp ; then
				taxid=$( grep "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp | head -n1 | cut -f1 )
				echo "taxid: "$taxid >> ${nucl_fa%.fasta}.tax_guide.blastx.out
				efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\n" -element TaxId,ScientificName,Rank >> ${nucl_fa%.fasta}.tax_guide.blastx.out
				sleep 0.4s
			fi
		fi
		if [ ! -s ${nucl_fa%.fasta}.tax_guide.blastx.out ] ; then
			echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out
		fi
		if grep -i -q "Caudovir\|Ackermannvir\|Herellevir\|Corticovir\|Levivir\|Tectivir\|crAss-like virus\|CrAssphage\|crassvirales\|Cyanophage\|Microvir\microphage\|Siphoviridae\|Myoviridae\|phage\|Podovir\|Halovir\|sphaerolipovir\|pleolipovir\|plasmid\|Inovir\|Ampullavir\|Bicaudavir\|Fusellovir\|Guttavir\|Ligamenvir\|Plasmavir\|Salterprovir\|Cystovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo ${nucl_fa%.fasta}.rotate.fasta >> DTR_seqs_for_phanotate.txt
		else
			echo ${nucl_fa%.fasta}.rotate.fasta >> DTR_seqs_for_prodigal.txt
		fi
	done
	if [ -s DTR_seqs_for_phanotate.txt ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running PHANOTATE, annotate DTR contigs " $MDYT
		cat DTR_seqs_for_phanotate.txt | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU phanotate.py -f fasta -o {}.rotate.phan.fasta {}.rotate.fasta
		for PHAN in *.phan.fasta ; do 
			if [ "$ENFORCE_START_CODON" == "True" ] ; then
				sed 's/ /@/g' ${PHAN} | bioawk -c fastx '{ print }' | awk '{ if ($2 ~ /^[ATCG]TG/) { print ">"$1 ; print $2 }}' | sed 's/@/ /g' > ${PHAN%.fasta}.sort.fasta
			else
				sed 's/ /@/g' ${PHAN} | bioawk -c fastx '{ print }' | awk '{ print ">"$1 ; print $2 }' | sed 's/@/ /g' > ${PHAN%.fasta}.sort.fasta
			fi
		done

		PHANQ=$( find . -maxdepth 1 -type f -regextype sed -regex ".*phan.sort.fasta" )
		echo "$PHANQ" | sed 's/\.phan\.sort\.fasta//g' | xargs -n 1 -I {} -P $CPU seqkit translate -x -T 11 {}.phan.sort.fasta -o {}.trans.fasta >/dev/null 2>&1
		
		for PHAN in *.phan.fasta ; do
			ORIG_CONTIG=$( grep ">" ${PHAN%.rotate.phan.fasta}.fna | cut -d " " -f2 ) ;
			sed 's/\[START=//g ; s/\]//g ; s/ \[SCORE=.*//g' ${PHAN%.phan.fasta}.trans.fasta | bioawk -v OC="$ORIG_CONTIG" -c fastx '{ split($1, start, "[. ]") ; split(start[2], sq, "[_]") ; if (substr($seq, 1, 1) == "M" || $4 > 3) {FIVE=""} else {FIVE="5primeInc"} ; if (substr($seq, length($seq), 1) == "*") {THREE=""} else {THREE="3primeInc"}; print ">" start[1]"_"NR " ["$4" - "sq[1]"] "FIVE THREE" "OC ; print $seq }' > ${PHAN%.phan.fasta}.AA.fasta
		done			
	fi
	if [ -s DTR_seqs_for_prodigal.txt ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running prodigal, annotate DTR contigs " $MDYT
		cat DTR_seqs_for_prodigal.txt | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU prodigal -a {}.rotate.prodigal.fasta -i {}.rotate.fasta -p meta -q >/dev/null 2>&1
		for PROD in *rotate.prodigal.fasta ; do
			sed 's/ /@/g' ${PROD} | bioawk -c fastx '{print}' | while read LINE ; do 
				ORIENTATION=$( echo "$LINE" | cut -d "#" -f 4 | sed 's/@//g' ) ;
				if [[ "$ORIENTATION" == 1 ]] ; then
					START_BASE=$( echo "$LINE" | cut -d "#" -f 2 | sed 's/@//g' ) ; 
					END_BASE=$( echo "$LINE" | cut -d "#" -f 3 | sed 's/@//g' ) ; 
				else
					START_BASE=$( echo "$LINE" | cut -d "#" -f 3 | sed 's/@//g' ) ; 
					END_BASE=$( echo "$LINE" | cut -d "#" -f 2 | sed 's/@//g' ) ; 
				fi	
				ORF_NAME=$( echo "$LINE" | cut -d "#" -f 1 | sed 's/@//g; s/\./_/g' ) ; 
				ORIG_CONTIG=$( grep ">" ${PROD%.prodigal.fasta}.fasta | cut -d " " -f2 ) ;
				AA_SEQ=$( echo "$LINE" | cut -f2 ) ;
				if echo $AA_SEQ | grep -q "\*" ; then
					INC3=""
				else
					INC3="3primeInc"
				fi
				FAA=${AA_SEQ:0:1}
				if [ "$FAA" != "M" ] && [ $START_BASE -le 3 ] ; then
					INC5="5primeInc"
				else
					INC5=""
				fi
				echo ">"${ORF_NAME} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${PROD%.prodigal.fasta}.AA.fasta
		done
	fi
	ROTATE_AAs=$( find . -maxdepth 1 -type f -name "${run_title}*rotate.AA.fasta" | sed 's/\.\///g' )
	if [ -n "$ROTATE_AAs" ] ; then
		for ROT in $ROTATE_AAs ; do 
			bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' $ROT > ${ROT%.fasta}.sorted.fasta
		done
	fi
fi

#-# 4 hhmscan circles/DTRs
ROTATE_SORT_AAs=$( find . -maxdepth 1 -type f -name "${run_title}*rotate.AA.sorted.fasta" | sed 's/\.\///g' )
if [ -n "$ROTATE_SORT_AAs" ] ; then
	#cat $( find . -maxdepth 1 -type f -name "${run_title}*rotate.AA.sorted.fasta" ) > all_DTR_sort_genome_proteins.AA.fasta
	for ROTQ in $ROTATE_SORT_AAs ; do
		cat $ROTQ
	done > all_DTR_sort_genome_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_DTR_sort_genome_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_sort_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_sort_genome_proteins.AA.fasta
		SPLIT_AA_DTR_sort=$( find . -maxdepth 1 -type f -name "SPLIT_DTR_sort_GENOME_AA_*.fasta" )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running hmmscan1, annotate DTR contigs " $MDYT
		if  [[ $virus_domain_db = "standard" ]] ; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
		elif [[ $virus_domain_db = "virion" ]]; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
		elif [[ $virus_domain_db = "rna_virus" ]]; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
		else
			echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
			rm -f ./*{0..9}.fasta
			break
		fi
		HMM_REP_NUMEBR=$( find . -maxdepth 1 -type f -name "SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan_replicate.out" | wc -l )
		if [[ $FOR_PLASMIDS = "True" ]]; then
			if [ $HMM_REP_NUMEBR -gt 0 ] ; then
				cat SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan.out SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_DTR_COMBINED.AA.hmmscan.sort.out
			else
				cat SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_DTR_COMBINED.AA.hmmscan.sort.out
			fi
		else
			if [ $HMM_REP_NUMEBR -gt 0 ] ; then
				cat SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan.out SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_DTR_COMBINED.AA.hmmscan.sort.out
			else
				cat SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_DTR_COMBINED.AA.hmmscan.sort.out
			fi
		fi		
		if [ -e SPLIT_DTR_COMBINED.AA.hmmscan.sort.out ] ; then
			cut -f3 SPLIT_DTR_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
				grep "${HIT}_" SPLIT_DTR_COMBINED.AA.hmmscan.sort.out > ${HIT}.rotate.AA.hmmscan.sort.out
				grep "${HIT}_" SPLIT_DTR_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.rotate.AA.called_hmmscan1.txt
				grep -v -f ${HIT}.rotate.AA.called_hmmscan1.txt ${HIT}.rotate.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.rotate.AA.no_hmmscan1.fasta
			done

		fi
	else
		echo "no AA seqs found for hmmscan1, annotate DTR contigs "
	fi
	for ROT in $ROTATE_SORT_AAs ; do 
		if [ ! -s ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan1.fasta ] ; then
			echo "all ORFs called for ${ROT%.rotate.AA.sorted.fasta} in first Hmmscan"
			#cp ${ROT} ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan1.fasta
		fi
	done	
	DTR_AA_FOR_HMM2=$( find . -maxdepth 1 -type f -name "${run_title}*AA.no_hmmscan1.fasta" )
	if [ -n "$DTR_AA_FOR_HMM2" ] ; then
		#cat $( find . -maxdepth 1 -type f -name "${run_title}*AA.no_hmmscan1.fasta" ) > all_DTR_HMM2_proteins.AA.fasta
		for DTRAQ in $DTR_AA_FOR_HMM2 ; do
			cat $DTRAQ
		done > all_DTR_HMM2_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_DTR_HMM2_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_HMM2_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_HMM2_proteins.AA.fasta
			SPLIT_DTR_HMM2=$( find . -maxdepth 1 -type f -name "SPLIT_DTR_HMM2_GENOME_AA_*.fasta" )
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running hmmscan2, annotate DTR contigs " $MDYT
			echo "$SPLIT_DTR_HMM2" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan2.out --cpu 1 -E 1e-8 --noali ${CENOTE_DBS}/hmmscan_DBs/useful_hmms_baits_and_not2a {}.fasta >/dev/null 2>&1
			cat SPLIT_DTR_HMM2_GENOME_AA_*AA.hmmscan2.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_DTR_HMM2_COMBINED.AA.hmmscan2.sort.out
		else
			echo "no AA seqs found for hmmscan2, annotate DTR contigs"
		fi
	fi
	if [ -s SPLIT_DTR_HMM2_COMBINED.AA.hmmscan2.sort.out ] ; then
		cut -f3 SPLIT_DTR_HMM2_COMBINED.AA.hmmscan2.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
			grep "${HIT}_" SPLIT_DTR_HMM2_COMBINED.AA.hmmscan2.sort.out > ${HIT}.rotate.AA.hmmscan2.sort.out
			grep "${HIT}_" SPLIT_DTR_HMM2_COMBINED.AA.hmmscan2.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.rotate.AA.called_hmmscan2.txt
			grep -v -f ${HIT}.rotate.AA.called_hmmscan2.txt ${HIT}.rotate.AA.no_hmmscan1.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.rotate.AA.no_hmmscan2.fasta
		done
	fi
#	for ROT in $ROTATE_SORT_AAs ; do 
#		if [ ! -s ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan2.fasta ] ; then
#			cp ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan1.fasta ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan2.fasta
#		fi
#	done
	for ROT_AAs in $ROTATE_SORT_AAs ; do
		echo ">Feature "${ROT_AAs%.rotate.AA.sorted.fasta}" Table1" > ${ROT_AAs%.rotate.AA.sorted.fasta}.SCAN.tbl
		if [ -s ${ROT_AAs%.AA.sorted.fasta}.AA.called_hmmscan1.txt ] || [ -s ${ROT_AAs%.AA.sorted.fasta}.AA.called_hmmscan2.txt ] ; then
			CALL_ALL_HMM=$( find . -maxdepth 1 -type f -regextype sed -regex "./${ROT_AAs%.rotate.AA.sorted.fasta}\..*called_hmmscan.*txt" | sed 's/\.\///g' )
			if [ -n "$CALL_ALL_HMM" ] ; then
				cat $( find . -maxdepth 1 -type f -regextype sed -regex "./${ROT_AAs%.rotate.AA.sorted.fasta}\..*called_hmmscan.*txt" | sed 's/\.\///g' ) > ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt

				if [ -s ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt ] ; then
					cat ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt | sed 's/ $//g' | while read LINE ; do 
						PROTEIN_INFO=$( grep "$LINE \[" ${ROT_AAs} ) ;  
						START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
						END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
						if [ -s ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan.sort.out ] && grep -q "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan.sort.out ; then
							HMM_INFO=$( grep "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' ) ; 
						else
							HMM_INFO=$( grep "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan2.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' )
						fi
						INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
						PROTEIN_NAME=$( echo $HMM_INFO | cut -d " " -f2- ) ; 
						echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""protein motif:$INFERENCEH" >> ${ROT_AAs%.rotate.AA.sorted.fasta}.SCAN.tbl ;
					done
				fi
			else
				echo "${ROT_AAs%.rotate.AA.sorted.fasta} (DTR contig) has no hits in hmmscan1 or hmmscan2 databases"
			fi
		fi
	done
fi


#- blastn all dtr seqs
if [ -n "$ROTATED_DTR_CONTIGS" ] && [ $handle_knowns == "blast_knowns" ] ; then
	if [ -s ${BLASTN_DB}.nsq ] || [ -s ${BLASTN_DB}.1.nsq ] || [ -s ${BLASTN_DB}.01.nsq ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running BLASTN, DTR contigs " $MDYT		


		### new blastn
		echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU blastn -task megablast -query {}.rotate.fasta -db ${BLASTN_DB} -outfmt '6 std qlen slen' -max_target_seqs 100 -perc_identity 90 -num_threads 1 -word_size 26 -evalue 1e-20 -out {}.blastn.out >/dev/null 2>&1
		for circle in $ROTATED_DTR_CONTIGS ; do
			if [ -s "${circle%.rotate.fasta}.blastn.out" ]; then
				python ${CENOTE_SCRIPT_DIR}/anicalc/anicalc.py -i ${circle%.rotate.fasta}.blastn.out -o ${circle%.rotate.fasta}.blastn_anicalc.out
				awk '{OFS="\t"}{FS="\t"}{ if (NR==1) {print $1, $2, $4, $5} else if ($4>=95 && $5>=85) {print $1, $2, $4, $5}}' ${circle%.rotate.fasta}.blastn_anicalc.out | head -n2 > ${circle%.rotate.fasta}.blastn_intraspecific.out
			fi
			if [ -s "${circle%.rotate.fasta}.blastn_intraspecific.out" ]; then
				INTRA_LINES=$( cat ${circle%.rotate.fasta}.blastn_intraspecific.out | wc -l | bc )	
				if [ "$INTRA_LINES" -ge 2 ] ; then
					#-#-#
					ORGANISM_H=$( head -n2 ${circle%.rotate.fasta}.blastn_intraspecific.out | tail -n1 | sed 's/\[/&\n/;s/.*\n//;s/\]/\n&/;s/\n.*//' )
					if grep -q "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp ; then
						taxid=$( grep "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp | head -n1 | cut -f1 )
						echo "taxid: "$taxid >> ${circle%.rotate.fasta}.tax_guide.blastn.out
						efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\n" -element TaxId,ScientificName,Rank >> ${circle%.rotate.fasta}.tax_guide.blastn.out
						sleep 0.4s
						efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -tab "\n" -element ScientificName >> ${circle%.rotate.fasta}.tax_guide.blastn.out
						sleep 0.4s
					fi

					if [ !  -z "${circle%.rotate.fasta}.tax_guide.blastn.out" ] ; then
						if grep -i -q "virus\|viridae\|virales\|Circular-genetic-element\|Circular genetic element\|plasmid\|phage" ${circle%.rotate.fasta}.tax_guide.blastn.out ; then
							#echo $circle " is closely related to a virus that has already been deposited in GenBank nt. "
							#cat ${circle%.rotate.fasta}.tax_guide.blastn.out
							cp ${circle%.rotate.fasta}.tax_guide.blastn.out ${circle%.rotate.fasta}.tax_guide.KNOWN_VIRUS.out

						else 
							#echo $circle "$(tput setaf 4) is closely related to a chromosomal sequence that has already been deposited in GenBank nt and will be checked for viral and plasmid domains. $(tput sgr 0)"
							#cat ${circle%.rotate.fasta}.tax_guide.blastn.out
							cp ${circle%.rotate.fasta}.tax_guide.blastn.out ${circle%.rotate.fasta}.tax_guide.CELLULAR.out
						fi
					fi
				fi
			fi
		done
		###

	else
		echo "BLASTN databases not found, skipping BLASTN step"
	fi
fi

PROTEIN_NO_HMMSCAN2=$( find . -maxdepth 1 -type f -name "*.rotate.AA.no_hmmscan2.fasta" )

if [ -n "$PROTEIN_NO_HMMSCAN2" ]; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running RPSBLAST, annotate DTR contigs" $MDYT
	#cat $( find . -maxdepth 1 -type f -name "*.rotate.AA.no_hmmscan2.fasta" ) > all_DTR_rps_proteins.AA.fasta
	for PROT_NOQ in $PROTEIN_NO_HMMSCAN2 ; do
		cat $PROT_NOQ
	done > all_DTR_rps_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_DTR_rps_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
	AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_RPS_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_rps_proteins.AA.fasta
		SPLIT_DTR_AA_RPS=$( find . -maxdepth 1 -type f -name "SPLIT_DTR_RPS_AA_*.fasta" )

		echo "$SPLIT_DTR_AA_RPS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t rpsblast -evalue 1e-4 -num_descriptions 5 -num_alignments 1 -db ${CENOTE_DBS}/cdd_rps_db/Cdd -seg yes -query {}.fasta -line_length 200 -out {}.rpsb.out >/dev/null 2>&1
		cat *rpsb.out | awk '{ if ($0 ~ /^>/) {printf $0 ; getline; print $0} else { print $0}}' > COMBINED_RESULTS.rotate.AA.rpsblast.out
		perl ${CENOTE_SCRIPT_DIR}/rpsblastreport2tbl_mt_annotation_pipe_biowulf.pl ;
		grep "protein_id	" COMBINED_RESULTS.NT.tbl | sed 's/.*protein_id	lcl|//g' | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read CONTIG ; do
			#echo "${CONTIG}"
			echo ">Feature ${CONTIG} Table1" > ${CONTIG}.NT.tbl

			grep -A2 -B1 "${CONTIG}_" COMBINED_RESULTS.NT.tbl | sed '/--/d' >> ${CONTIG}.NT.tbl
		done
	else
		echo "no AA seqs for RPSBLAST, annotate DTR contigs"
	fi
fi

# Detecting any tRNAs and making a tbl addenum file
if [ -n "$ROTATED_DTR_CONTIGS" ]; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running tRNAscan-SE " $MDYT
	for GENOME_NAME in $ROTATED_DTR_CONTIGS ; do
		tRNAscan-SE -Q -G -o ${GENOME_NAME}.trnascan-se2.txt ${GENOME_NAME} >/dev/null 2>&1
		
		if grep -q "${GENOME_NAME%.rotate.fasta}" ${GENOME_NAME}.trnascan-se2.txt ;then
			#echo "$(tput setaf 5) "$GENOME_NAME" was found to encode tRNA(s); making .tbl file $(tput sgr 0)"

			grep "${GENOME_NAME%.rotate.fasta}" $GENOME_NAME.trnascan-se2.txt | while read LINE ; do 
				TRNA_START=$( echo $LINE | cut -d " " -f3 ) ; 
				TRNA_END=$( echo $LINE | cut -d " " -f4 ) ; 
				TRNA_NUMBER=$( echo $LINE | cut -d " " -f2 ) ; 
				TRNA_TYPE=$( echo $LINE | cut -d " " -f5 ) ; 
				TRNA_SCORE=$( echo $LINE | cut -d " " -f9 ) ; 
				echo -e "$TRNA_START\t""$TRNA_END\t""tRNA\n""\t\t\tgene\t""${GENOME_NAME%.rotate.fasta}""_tRNA$TRNA_NUMBER\n""\t\t\tproduct\t""tRNA-$TRNA_TYPE\n""\t\t\tinference\t""tRNAscan-SE score:$TRNA_SCORE" >> ${GENOME_NAME%.rotate.fasta}.trna.tbl; 
			done
		fi
	done
fi

NEW_FASTAS=$( echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate//g' )

if [ -n "$NEW_FASTAS" ]; then
	for nucl_fa in $NEW_FASTAS ; do
		if [ -s ${nucl_fa%.fasta}.NT.tbl ] && [ -s ${nucl_fa%.fasta}.SCAN.tbl ] && [ -s ${nucl_fa%.fasta}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fasta}.BLASTP.tbl > ${nucl_fa%.fasta}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			grep -v -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' ${nucl_fa%.fasta}.NT.tbl | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' >> ${nucl_fa%.fasta}.int.tbl ;
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			cat ${nucl_fa%.fasta}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.SCAN.tbl ] && [ -s ${nucl_fa%.fasta}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fasta}.BLASTP.tbl > ${nucl_fa%.fasta}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			cat ${nucl_fa%.fasta}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.NT.tbl ] && [ -s ${nucl_fa%.fasta}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fasta}.BLASTP.tbl > ${nucl_fa%.fasta}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			grep -v -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' ${nucl_fa%.fasta}.NT.tbl | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' >> ${nucl_fa%.fasta}.int.tbl ;
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.NT.tbl ] && [ -s ${nucl_fa%.fasta}.SCAN.tbl ] ; then
			cat ${nucl_fa%.fasta}.NT.tbl > ${nucl_fa%.fasta}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			cat ${nucl_fa%.fasta}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.NT.tbl ] ; then
			cat ${nucl_fa%.fasta}.NT.tbl > ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.SCAN.tbl ] ; then

			cat ${nucl_fa%.fasta}.SCAN.tbl >> ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fasta}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fasta}.BLASTP.tbl > ${nucl_fa%.fasta}.int.tbl
			if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
				cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
			fi
		else
			echo ">Feature ${nucl_fa%.fasta} Table1" > ${nucl_fa%.fasta}.int.tbl
		fi
	done
fi

# remove ORFs within ORFs that are 'hypothetical'
if [ "$ORF_WITHIN" == "True" ] ; then
	echo "$(tput setaf 5) Removing ORFS within ORFs that are 'hypothetical' $(tput sgr 0)"

	INT_TBL=$( find . -maxdepth 1 -type f -name "*.int.tbl" )
	if [ -n "$INT_TBL" ] ; then
		for feat_tbl3 in $INT_TBL ; do
			grep "^[0-9]" $feat_tbl3 | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.int.tbl}.all_start_stop.txt ;
			cat "${feat_tbl3%.int.tbl}.all_start_stop.txt" | while read linev ; do
				all_start=$( echo $linev | cut -d " " -f1 )
				all_end=$( echo $linev | cut -d " " -f2 )
				if [[ "$all_end" -gt "$all_start" ]]; then
					for ((counter_f=(( $all_start + 1 ));counter_f<=$all_end;counter_f++)); do
						echo " " "$counter_f" " " >> ${feat_tbl3%.int.tbl}.used_positions.txt
						
					done
				elif [[ "$all_start" -gt "$all_end" ]]; then
					for ((counter_r=$all_end;counter_r<=(( $all_start - 1 ));counter_r++)) ; do
						echo " " "$counter_r" " " >> ${feat_tbl3%.int.tbl}.used_positions.txt
					done
				fi
			done
			sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g ; s/length=.*//g' $feat_tbl3 | sed '/--/d' > ${feat_tbl3%.int.tbl}.comb2.tbl ; 
			grep -e 'hypothetical protein' -B2 ${feat_tbl3%.int.tbl}.comb2.tbl | grep "^[0-9]" | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.int.tbl}.hypo_start_stop.txt ;

			# Remove redudant ORFs that are subORFs of ORFs overlapping the wrap-point
			GENOME_LENGTH=$( bioawk -c fastx '{print length($seq)}' ${feat_tbl3%.int.tbl}.rotate.fasta )
			cat "${feat_tbl3%.int.tbl}.all_start_stop.txt" | while read linet ; do
				all_start=$( echo $linet | cut -d " " -f1 )
				all_end=$( echo $linet | cut -d " " -f2 )
				if [[ "$all_end" -gt "$GENOME_LENGTH" ]] ; then
					cat "${feat_tbl3%.int.tbl}.all_start_stop.txt" | while read lineq ; do
						q_end=$( echo $lineq | cut -d " " -f2 )
						OVERWRAP_END=$(( $q_end + $GENOME_LENGTH ))
						if [[ "$all_end" = "$OVERWRAP_END"  ]] ; then
							echo "$lineq" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
						fi
					done
				elif [[ "$all_start" -gt "$GENOME_LENGTH" ]] ; then
					cat "${feat_tbl3%.int.tbl}.all_start_stop.txt" | while read lineq ; do
						q_start=$( echo $lineq | cut -d " " -f2 )
						OVERWRAP_START=$(( $q_start + $GENOME_LENGTH ))
						if [[ "$all_start" = "$OVERWRAP_START"  ]] ; then
							echo "$lineq" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
						fi
					done
				fi
			done


			cat "${feat_tbl3%.int.tbl}.hypo_start_stop.txt" | while read liney ; do
				loc_start=$( echo $liney | cut -d " " -f1 )
				loc_end=$( echo $liney | cut -d " " -f2 )
				loc1_start=$( echo " " "$loc_start" " " )
				if grep -q "$loc1_start" ${feat_tbl3%.int.tbl}.used_positions.txt ; then 

					if [[ "$loc_end" -gt "$loc_start" ]]; then
						gen_len=$(( $loc_end - $loc_start ))

						if [[ "$gen_len" -gt 1000 ]]; then
							continue
						else
							f_end=$(( $loc_end + 1 ))
							f1_end=$( echo " " "$f_end" " ")
							if grep -q "$f1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then

								echo "$liney" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
							fi
						fi
					else
						gen_len=$(( $loc_start - $loc_end ))

						if [[ "$gen_len" -gt 1000 ]]; then
							continue
						else
						r_end=$(( $loc_end - 1 ))
						r1_end=$( echo " " "$r_end" " ")

							if grep -q "$r1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then

									echo "$liney" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
							fi
						fi
					fi
				fi
			done
			if [ -s "${feat_tbl3%.int.tbl}.remove_hypo.txt" ]; then
				grep ">Feature" ${feat_tbl3%.int.tbl}.comb2.tbl | sed '/--/d' > ${feat_tbl3%.int.tbl}.int2.tbl
				grep -v -f ${feat_tbl3%.int.tbl}.remove_hypo.txt ${feat_tbl3%.int.tbl}.comb2.tbl | grep "^[0-9]" -A3 | sed '/--/d' >> ${feat_tbl3%.int.tbl}.int2.tbl ;
				else
					cp ${feat_tbl3%.int.tbl}.comb2.tbl ${feat_tbl3%.int.tbl}.int2.tbl
			fi
		done
	fi
else
	INT_TBL=$( find . -maxdepth 1 -type f -name "*.int.tbl" )
	if [ -n "$INT_TBL" ] ; then
		for feat_tbl3 in $INT_TBL ; do
			sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g ; s/length=.*//g' $feat_tbl3 | sed '/--/d' > ${feat_tbl3%.int.tbl}.int2.tbl ; 
		done
	fi
fi

# Grabbing ORFs wihout RPSBLAST hits and separating them into individual files for HHsearch
echo "$(tput setaf 5) Grabbing ORFs wihout RPS-BLAST hits and separating them into individual files for HHsearch $(tput sgr 0)"

INT2_TBL=$( find . -maxdepth 1 -type f -name "*.int2.tbl" )
if [ -n "$INT2_TBL" ] ; then
	for blastp_tbl1 in $INT2_TBL ; do
		grep -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' -B2 $blastp_tbl1 | grep "^[0-9]" | awk '{print $1 " - " $2}' > ${blastp_tbl1%.int2.tbl}.for_hhpred.txt ;
		grep -f ${blastp_tbl1%.int2.tbl}.for_hhpred.txt -A1 ${blastp_tbl1%.int2.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${blastp_tbl1%.int2.tbl}.rotate.blast_hypo.fasta ;
		csplit -z ${blastp_tbl1%.int2.tbl}.rotate.blast_hypo.fasta '/>/' '{*}' --prefix=${blastp_tbl1%.int2.tbl}.rotate --suffix-format=%02d.for_hhpred.fasta >/dev/null 2>&1
	done
fi

# Running HHsearch on remaining ORFs
if  [[ $HHSUITE_TOOL = "hhsearch" ]] || [[ $HHSUITE_TOOL = "hhblits" ]] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running HHsearch or HHblits " $MDYT
fi
dark_orf_list=$( find . -maxdepth 1 -type f -name "*.for_hhpred.fasta" )
#-#- can this be parallelized?
if [ -n "$dark_orf_list" ] ; then
	if  [[ $HHSUITE_TOOL = "hhsearch" ]] ; then
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU hhsearch -i {}.for_hhpred.fasta "${HHSUITE_DB_STR}" -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	elif [[ $HHSUITE_TOOL = "hhblits" ]] ; then
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU hhblits -i {}.for_hhpred.fasta "${HHSUITE_DB_STR}" -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	else
		echo "$(tput setaf 5) Valid option for HHsuite tool (i.e. hhsearch or hhblits) was not provided. Skipping step for "$dark_orf" $(tput sgr 0)"
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	fi
	cat_list=$( find . -maxdepth 1 -type f -name "*.out.hhr" )
	if [ -n "$cat_list" ] ; then
		cat *out.hhr > ${run_title}.rotate.out_all.hhr
		rm -f *out.hhr
	fi	
fi

rm -f *.rotate.AA.fasta

perl ${CENOTE_SCRIPT_DIR}/hhpredreport2tbl_mt_annotation_pipe_biowulf1_gjs_edits.pl

HH_TBL=$( find . -maxdepth 1 -type f -name "*.HH.tbl" )
if [ -n "$HH_TBL" ] ; then
	for HH_tbl1 in $HH_TBL ; do 
		sed 's/OS=.*//g; s/ ;//g; s/similar to AA sequence:UniProtKB:>\([0-9][A-Z].*\)/protein motif:PDB:\1/g; s/UniProtKB:>tr|.*|\(.\)/UniProtKB:\1/g; s/similar to AA sequence:UniProtKB:>\([a-z].*\)/protein motif:Scop:\1/g; s/similar to AA sequence:UniProtKB:>\(PF.*\)/protein motif:PFAM:\1/g; s/ is .*//g; s/ are .*//g' $HH_tbl1 | sed '/product/ s/; [a-zA-Z0-9_]\{1,20\}//g; s/;.*//g' > ${HH_tbl1%.HH.tbl}.HH2.tbl
	done
fi

# Combining tbl files from all search results AND fix overlapping ORF module
echo "$(tput setaf 5) Combining tbl files from all search results AND fix overlapping ORF module $(tput sgr 0)"

if [ -n "$INT2_TBL" ] ; then
	for feat_tbl4 in $INT2_TBL ; do 
		if [ -s "${feat_tbl4%.int2.tbl}.HH2.tbl" ] && [ -s "$feat_tbl4" ] ; then
			head -n1 $feat_tbl4 > ${feat_tbl4%.int2.tbl}.comb3.tbl
			grep -v -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' -e 'putative phage protein' $feat_tbl4 | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d; s/TPA_asm: //g; s/TPA://g' >> ${feat_tbl4%.int2.tbl}.comb3.tbl
			sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g; s/[bB]rain cDNA.*/hypothetical protein/g; s/(E\.C.*//g ; s/length=.*//g' ${feat_tbl4%.int2.tbl}.HH2.tbl | grep -v ">Feature" | sed '/--/d' >> ${feat_tbl4%.int2.tbl}.comb3.tbl ;
			if [ -s "${feat_tbl4%.int2.tbl}.ITR.tbl" ] ; then
				cat ${feat_tbl4%.int2.tbl}.ITR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
			fi
			if [ -s "${feat_tbl4%.int2.tbl}.DTR.tbl" ] ; then
				cat ${feat_tbl4%.int2.tbl}.DTR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
			fi 
		else
			cat $feat_tbl4 | sed '/--/d; s/TPA_asm: //g; s/TPA://g' > ${feat_tbl4%.int2.tbl}.comb3.tbl
			if [ -s "${feat_tbl4%.int2.tbl}.ITR.tbl" ] ; then
				cat ${feat_tbl4%.int2.tbl}.ITR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
			fi
			if [ -s "${feat_tbl4%.int2.tbl}.DTR.tbl" ] ; then
				cat ${feat_tbl4%.int2.tbl}.DTR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
			fi 
		fi
		GENOME_LENGTH=$( bioawk -c fastx '{print length($seq)}' ${feat_tbl4%.int2.tbl}.rotate.fasta )
		grep "^[0-9]" ${feat_tbl4%.int2.tbl}.comb3.tbl | cut -f1,2 | while read LINE ; do 
			START_SITE=$( echo $LINE | awk -v TOTAL_LENGTH="$GENOME_LENGTH" '{ if ($1>TOTAL_LENGTH) { print $1 } else {print "not_overlaps"}}' ) ; 
			END_SITE=$( echo $LINE | awk -v TOTAL_LENGTH="$GENOME_LENGTH" '{ if ($2>TOTAL_LENGTH) { print $2 } else {print "not_overlaps"}}' ) ;
			if [ "$START_SITE" != 'not_overlaps' ] || [ "$END_SITE" != 'not_overlaps' ] ; then
	#			echo $START_SITE
	#			echo $END_SITE
	#			echo $GENOME_LENGTH
				if [ "$START_SITE" != 'not_overlaps' ] ; then
					END_SITEQ=$( echo $LINE | awk '{ print $2 }' )				
					NEW_START=$(( $START_SITE - $GENOME_LENGTH ))
					FIRST_HALF=$( cat ${feat_tbl4%.int2.tbl}.comb3.tbl | sed "/$START_SITE	$END_SITEQ	CDS/q" | sed "s/$START_SITE	$END_SITEQ	CDS/$NEW_START	<1	CDS/g" )
					#echo "$FIRST_HALF"
					SECOND_HALF=$( cat ${feat_tbl4%.int2.tbl}.comb3.tbl | sed "s/$START_SITE	$END_SITE	CDS/$GENOME_LENGTH	$END_SITE/g" | sed -n '/</,$p' ) ;
				elif [ "$END_SITE" != 'not_overlaps' ] ; then
					START_SITEQ=$( echo $LINE | awk '{ print $1 }' )
					NEW_END=$(( $END_SITE - $GENOME_LENGTH ))
					FIRST_HALF=$( cat ${feat_tbl4%.int2.tbl}.comb3.tbl | sed "/$START_SITEQ	$END_SITE	CDS/q" | sed "s/$START_SITEQ	$END_SITE	CDS/$START_SITEQ	$GENOME_LENGTH	CDS/g" )
					SECOND_HALF=$( cat ${feat_tbl4%.int2.tbl}.comb3.tbl | sed "s/$START_SITEQ	$END_SITE	CDS/<1	$NEW_END/g" | sed -n '/</,$p' )
				fi
				echo -e "$FIRST_HALF\n""$SECOND_HALF" > tempq.${feat_tbl4%.int2.tbl}.comb3.tbl
				mv tempq.${feat_tbl4%.int2.tbl}.comb3.tbl ${feat_tbl4%.int2.tbl}.comb3.tbl

			fi

		done ; 

	done
fi
### insert comb3.tbl edits
COMB3_TBL=$( find . -maxdepth 1 -type f -name "*.comb3.tbl" )
if [ -n "$COMB3_TBL" ] ; then
	for comb3 in $COMB3_TBL ; do
		## tRNA overlap
		if grep -q "[0-9]	tRNA" $comb3 ; then
			grep "[0-9]	tRNA" $comb3 | awk -v name="${comb3%.comb3.tbl}" '{OFS="\t"}{FS="\t"}{ if ($1<$2) {print name, $1, $2, "fwd"} else {print name, $2, $1, "rev"}}' > ${comb3%.comb3.tbl}.tRNA.bed
		fi
		if grep -q "	CDS" $comb3 ; then
			grep "	CDS" $comb3 | awk -v name="${comb3%.comb3.tbl}" '{OFS="\t"}{FS="\t"}{ if ($1<$2) {print name, $1, $2, "fwd"} else {print name, $2, $1, "rev"}}' > ${comb3%.comb3.tbl}.CDS.bed
		fi
		if [ -s ${comb3%.comb3.tbl}.CDS.bed ] && [ -s ${comb3%.comb3.tbl}.tRNA.bed ] ; then
			bedtools intersect -wa -a ${comb3%.comb3.tbl}.CDS.bed -b ${comb3%.comb3.tbl}.tRNA.bed | awk '{OFS="\t"}{FS="\t"}{ if ($4=="fwd") {print $2, $3} else {print $3, $2}}' > ${comb3%.comb3.tbl}.ORFs_over_tRNAs.tsv
			if [ -s ${comb3%.comb3.tbl}.ORFs_over_tRNAs.tsv ] ; then
				head -n1 ${comb3} > ${comb3}.tmp
				grep -v -f ${comb3%.comb3.tbl}.ORFs_over_tRNAs.tsv ${comb3} | grep -A3 "^[0-9]" >> ${comb3}.tmp
				sed '/--/d' ${comb3}.tmp > ${comb3}
			fi
		fi

		##
		if grep -q "5primeInc\|3primeInc" ${comb3%.comb3.tbl}.rotate.AA.sorted.fasta ; then
			grep "5primeInc\|3primeInc" ${comb3%.comb3.tbl}.rotate.AA.sorted.fasta | while read INCOMPLETE ; do
				START_BASEH=$( echo $INCOMPLETE | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
				END_BASEH=$( echo $INCOMPLETE | sed 's/.*- \(.*\)\].*/\1/' ) ; 
				if echo "$INCOMPLETE" | grep -q "5primeInc" && echo "$INCOMPLETE" | grep -q "3primeInc" ; then
					if grep -q "^${START_BASEH}	${END_BASEH}" $comb3 ; then
						sed -i -e ':a' -e 'N' -e '$!ba' -e "s/${START_BASEH}	${END_BASEH}	CDS\n/<${START_BASEH}	>${END_BASEH}	CDS\n			codon_start	1\n/g" $comb3
					fi
				elif echo "$INCOMPLETE" | grep -q "5primeInc" ; then
					if grep -q "^${START_BASEH}	${END_BASEH}" $comb3 ; then
						sed -i -e ':a' -e 'N' -e '$!ba' -e "s/${START_BASEH}	${END_BASEH}	CDS\n/<${START_BASEH}	${END_BASEH}	CDS\n			codon_start	1\n/g" $comb3
					fi
				elif echo "$INCOMPLETE" | grep -q "3primeInc" ; then
					if grep -q "^${START_BASEH}	${END_BASEH}" $comb3 ; then
						sed -i -e ':a' -e 'N' -e '$!ba' -e "s/${START_BASEH}	${END_BASEH}	CDS\n/${START_BASEH}	>${END_BASEH}	CDS\n			codon_start	1\n/g" $comb3
					fi
				fi
			done
		fi
		## bad names fix
		sed -i 's/product	 /product	/g ; s/product	-/product	/g ; s/product	\;/product	/g ; s/product	=/product	/g ; s/product	_/product	/g ; s/product	\; /product	/g ; s/Length=.*//g ; s/\; Provisional.//g ; s/\; Validated.//g ; s/: .*//g ; s/\; Reviewed.//g ; s/\;$//g' $comb3
		##
	done
fi

###
rm -f *tmp.tbl
COMB3_TBL=$( find . -maxdepth 1 -type f -name "*.comb3.tbl" )
if [ -n "$COMB3_TBL" ] ; then
	for feat_tbl2 in *.comb3.tbl ; do 
		if grep -i -q "CRESS\|genomovir\|circovir\|bacilladnavir\|redondovir\|nanovir\|geminivir\|smacovir" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
			echo ${feat_tbl2%.comb3.tbl}" is a CRESS virus of some kind"
		else
			CONJ_COUNT=$( grep -i "virb\|type-IV\|secretion system\|conjugation\|conjugal\|transposon\| tra[a-z] \|	tra[a-z]\|trb[b-z]\|pilus" $feat_tbl2 | grep -v "TRAF\|TRAP\|protein tyrosine phosphatase\|ttRBP\|SpoU\|transport\|central" | wc -l )
			STRUCTURAL_COUNT=$( grep -i "capsid\|terminase\|portal\|baseplate\|base plate\|tail\|collar\|zot\|zonular\|minor coat\|packaging\|	virion protein" $feat_tbl2 | wc -l )
			if [[ $CONJ_COUNT -gt 0 ]] && [[ $STRUCTURAL_COUNT == 0 ]] ; then
				TAX_ORF="Conjugative Transposon"
			elif grep -i -q "large terminase\|large subunit terminase\|packaging\|terminase, large\|terminase large" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "large terminase\|large subunit terminase\|packaging\|terminase, large\|terminase large" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )
			elif grep -i -q "dnab\|dna polymerase\|polb\|rdrp\|rna dependent rna polymerase" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "dnab\|dna polymerase\|polb\|rdrp\|rna dependent rna polymerase" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )			
			elif grep -i -q "portal" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "portal" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )		
			elif grep -i -q "rep \|replica\|repa " $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "rep \|replica\|repa" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )	
			elif grep -i -q "capsid\|cap \|stnv\|coat" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "capsid\|cap \|stnv\|coat" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )
			elif grep -i -q "baseplate\|base-plate\|base plate" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "baseplate\|base-plate\|base plate" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )
			elif grep -i -q "tail" $feat_tbl2 ; then
				TAX_ORF=$( grep -i -B1 "tail" $feat_tbl2 | head -n1 | sed 's/.*lcl|\(.*\)/\1/' )
			else 
				TAX_ORF="No_suitable_orf"
			fi
			if [ "$TAX_ORF" == "Conjugative Transposon" ] ; then
				#echo "${feat_tbl2%.comb3.tbl} looks like a conjugative transposon"
				echo $TAX_ORF > ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out

			elif [ "$TAX_ORF" == "No_suitable_orf" ] ; then
				echo "No suitable ORF for taxonomy found for ${feat_tbl2%.comb3.tbl}, using BLASTX result."
			else

				grep -A1 "$TAX_ORF " ${feat_tbl2%.comb3.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${feat_tbl2%.comb3.tbl}.tax_orf.fasta
				if [[ $STRUCTURAL_COUNT == 1 ]] || [[ $STRUCTURAL_COUNT -gt 1 ]] ; then
					blastp -evalue 1e-2 -outfmt "6 qseqid stitle pident evalue length" -num_threads $CPU -num_alignments 1 -db ${CENOTE_DBS}/blast_DBs/virus_adinto_polinton_prot_190925 -query ${feat_tbl2%.comb3.tbl}.tax_orf.fasta -out ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out >/dev/null 2>&1
				else
					blastp -evalue 1e-2 -outfmt "6 qseqid stitle pident evalue length" -num_threads $CPU -num_alignments 1 -db ${CENOTE_DBS}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query ${feat_tbl2%.comb3.tbl}.tax_orf.fasta -out ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out >/dev/null 2>&1
				fi
				if [ ! -s "${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out" ]; then
					echo "unclassified virus" > ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ;
				elif grep -q "virophage" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Virophage	Unclassified Taxon" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
				elif grep -q "adinto" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Adintovirus	Unclassified Taxon" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
				elif grep -i -q "polinton" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Polinton-like virus	Unclassified Taxon" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out				
				else
					#-#-#
					ORGANISM_H=$( head -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out | sed 's/\[/&\n/;s/.*\n//;s/\]/\n&/;s/\n.*//' )
					if grep -q "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp ; then
						taxid=$( grep "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp | head -n1 | cut -f1 )
						echo "taxid: "$taxid >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
						efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\n" -element TaxId,ScientificName,Rank >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
						sleep 0.4s
					fi

				fi
			fi
		fi
	done
fi

# module for taxonomy of all hallmark genes
cd ${base_directory}/${run_title}

if [ "$HALLMARK_TAX" == "True" ] && [ -d DTR_contigs_with_viral_domain/ ] ;then
	HALLMARK_FILES=$( find DTR_contigs_with_viral_domain/ -maxdepth 1 -type f -name "*rotate.AA.hmmscan.sort.out" | sed 's/\.\///g' )
	if [ -n "${HALLMARK_FILES}" ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: reporting taxonomy for each hallmark gene, circular contigs " $MDYT
		echo "$HALLMARK_FILES" | while read SCAN ; do
			cut -f3 $SCAN | while read HALLMARK ; do 
				seqkit grep -p "$HALLMARK" ${SCAN%.hmmscan.sort.out}.sorted.fasta ; 
			done > ${SCAN%.hmmscan.sort.out}.hallmarks.fasta

			blastp -evalue 1e-2 -outfmt "6 qseqid stitle pident evalue length" -num_threads $CPU -num_alignments 1 -db ${CENOTE_DBS}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query ${SCAN%.hmmscan.sort.out}.hallmarks.fasta -out ${SCAN%.hmmscan.sort.out}.hallmarks.blastp.out >/dev/null 2>&1

			sort -k1,1 ${SCAN%.hmmscan.sort.out}.hallmarks.blastp.out | while read LINE ; do 
				ORGANISM_H=$( echo "$LINE" | sed 's/\[/&\n/;s/.*\n//;s/\]/\n&/;s/\n.*//' )
				if grep -q "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp ; then
					echo $LINE 
					echo "taxid: "$taxid
					taxid=$( grep "	|	${ORGANISM_H}	|	" ${CENOTE_DBS}/taxdump/names.dmp | head -n1 | cut -f1 )
					efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -block "*/Taxon" -tab "\n" -element TaxId,ScientificName,Rank 
					echo "@@@@@@"
					sleep 0.4s
				fi
			done > ${SCAN%.hmmscan.sort.out}.hallmarks.taxonomy.out
		done
	fi
fi


LIST_OF_ITR_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "ITR_containing_contigs/*fna" )

if [ -n "$LIST_OF_ITR_DOMAIN_CONTIGS" ] ; then
	if [ ! -d no_end_contigs_with_viral_domain ] ; then
		mkdir no_end_contigs_with_viral_domain
	fi
	cd ITR_containing_contigs
	LIST_OF_ITR_DOMAIN_CONTIGS=$( find . -maxdepth 1 -type f -regextype sed -regex ".*.fna" )
	echo "$LIST_OF_ITR_DOMAIN_CONTIGS" | sed 's/.fna//g ; s/\.\///g' | while read ITR_SEQ ; do
		if [ "$PROPHAGE" == "True" ] ;then
			sed "s/${ITR_SEQ}/${ITR_SEQ}_vs99/g" ${ITR_SEQ}.fna > ../no_end_contigs_with_viral_domain/${ITR_SEQ}_vs99.fna
			grep -A1 "^[0-9]\{1,6\}	[0-9]\{1,6\}	repeat_region" ${ITR_SEQ}.ITR.tbl | sed '/--/d' > ../no_end_contigs_with_viral_domain/${ITR_SEQ}_vs99.ITR.tbl			
		else
			cp ${ITR_SEQ}.fna ../no_end_contigs_with_viral_domain/${ITR_SEQ}.fna
			grep -A1 "^[0-9]\{1,6\}	[0-9]\{1,6\}	repeat_region" ${ITR_SEQ}.ITR.tbl | sed '/--/d' > ../no_end_contigs_with_viral_domain/${ITR_SEQ}.ITR.tbl			
		fi
	done
	cd ${base_directory}/${run_title}
else
	echo "No ITR contigs with minimum hallmark genes found."
fi

LIST_OF_VIRAL_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "no_end_contigs_with_viral_domain/*fna" )

if [ -n "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ; then
	. ${CENOTE_SCRIPT_DIR}/annotate_linear_contigs_2.1.5.sh
else
	echo "No linear contigs with minimum hallmark genes found."
fi

cd ${base_directory}/${run_title}

CIRCULAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -wholename "DTR_contigs_with_viral_domain/*fna" | grep -v "\.rc\.fna" )
if [ -n "$CIRCULAR_HALLMARK_CONTIGS" ] ; then
	for CIRC in $CIRCULAR_HALLMARK_CONTIGS ; do
		sed 's/ /#/g' $CIRC | bioawk -c fastx '{print ">"$name" DTR" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
	done
fi
if [ "$PROPHAGE" == "True" ] ;then
	LINEAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -regextype sed -regex ".*_vs[0-9]\{1,2\}.fna" )
	if [ -n "$LINEAR_HALLMARK_CONTIGS" ] ; then
		for LIN in $LINEAR_HALLMARK_CONTIGS ; do
			if [ -s ${LIN%.fna}.ITR.tbl ] ; then
				ORIGINAL_NAME=$( head -n1 ${LIN} | cut -d " " -f2 )
			else
				ORIGINAL_NAME=$( head -n1 ${LIN%_vs[0-9][0-9].fna}.fna | cut -d " " -f2 )
			fi
			if [ -s ${LIN%.fna}.ITR.tbl ] ; then
				sed 's/ /#/g' $LIN | bioawk -v ORI="$ORIGINAL_NAME" -c fastx '{print ">"$name" "ORI" ITR" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			else
				sed 's/ /#/g' $LIN | bioawk -v ORI="$ORIGINAL_NAME" -c fastx '{print ">"$name" "ORI" no_end_feature" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			fi
		done
	fi
else
	if [ -n "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ; then
		for LIN in $LIST_OF_VIRAL_DOMAIN_CONTIGS ; do
			if [ -s ${LIN%.fna}.ITR.tbl ] ; then
				sed 's/ /#/g' $LIN | bioawk -c fastx '{print ">"$name" ITR" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			else
				sed 's/ /#/g' $LIN | bioawk -c fastx '{print ">"$name" no_end_feature" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			fi
		done
	fi
fi		

cd ${base_directory}/${run_title}

COMB3_TBL=$( find * -maxdepth 1 -type f -name "*.comb3.tbl" )
if [ -n "$COMB3_TBL" ] ; then
	if [ ! -d "sequin_and_genome_maps" ]; then
		mkdir sequin_and_genome_maps
	fi
	for feat_tbl2 in $COMB3_TBL ; do
		JUST_TBL2_FILE=$( echo "$feat_tbl2" | sed 's/.*\///g' )
		file_core=${JUST_TBL2_FILE%.comb3.tbl}
		#echo $file_core
		file_numbers=$( echo ${file_core#${run_title}} | sed 's/[a-zA-Z]//g' )

		#echo $file_numbers
		tax_info=${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
		#echo $tax_info
		#-#-#
		TAX_EVALUE=$( head -n1 $tax_info | cut -f4 )
		TAX_PERCID=$( head -n1 $tax_info | cut -f3 )
		TAX_CUTOFF=$( echo -e "${TAX_EVALUE}\t${TAX_PERCID}" | awk '{if ($1<1e-100 && $2>90) {print "genus"} else if ($1<1e-20 && $2>40) {print "family"} else if ($1<1e-4 && $2>25) {print "order"} else {print "unclassified"}}' )
		if [ $TAX_CUTOFF == "genus" ] ; then
			if grep -q "Virophage	Unclassified Taxon" $tax_info ; then
				vir_name="Virophage" ;
			elif grep -q "Adintovirus	Unclassified Taxon" $tax_info ; then
				vir_name="Adintovirus" ;
			elif grep -q "Polinton-like virus	Unclassified Taxon" $tax_info ; then
				vir_name="Polinton-like virus" ;
			elif grep -q "	genus$" $tax_info ; then
				vir_name=$( grep "	genus$" $tax_info | cut -f2 )
			elif grep -q "	family$" $tax_info ; then
				vir_name=$( grep "	family$" $tax_info | cut -f2 )
			elif grep -q "	order$" $tax_info ; then
				vir_name=$( grep "	order$" $tax_info | cut -f2 )
			elif grep -q "Conjugative Transposon" $tax_info ; then
				vir_name="Conjugative Transposon" ;
			elif grep -q "No homologues found" $tax_info ; then
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="genetic element" ;
				else
					vir_name="circular genetic element" ;
				fi
			elif grep -q "Circular genetic element" $tax_info ; then
				vir_name="Circular genetic element" ;
			elif grep -q "dsRNA virus" $tax_info ; then
				vir_name="dsRNA virus" ;
			elif grep -q "ssRNA virus" $tax_info ; then
				vir_name="ssRNA virus" ;
			elif grep -q "unclassified RNA virus" $tax_info ; then
				vir_name="unclassified RNA virus" ;
			elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
				vir_name="unclassified ssDNA bacterial virus" ;
			elif grep -q -i "phage" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "plasmid" $tax_info ; then
				vir_name="metagenomic plasmid" ;
			elif grep -q "Bacteria" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "unclassified virus" $tax_info ; then
				vir_name="virus" ;		
			elif grep -q -i "virus\|viridae" $tax_info ; then
				vir_name="virus" ;
			else
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="unclassified element" ;
				else
					vir_name="Circular genetic element" ;
				fi
			fi		
		elif [ $TAX_CUTOFF == "family" ] ; then
			if grep -q "Virophage	Unclassified Taxon" $tax_info ; then
				vir_name="Virophage" ;
			elif grep -q "Adintovirus	Unclassified Taxon" $tax_info ; then
				vir_name="Adintovirus" ;
			elif grep -q "Polinton-like virus	Unclassified Taxon" $tax_info ; then
				vir_name="Polinton-like virus" ;
			elif grep -q "	family$" $tax_info ; then
				vir_name=$( grep "	family$" $tax_info | cut -f2 )
			elif grep -q "	order$" $tax_info ; then
				vir_name=$( grep "	order$" $tax_info | cut -f2 )
			elif grep -q "Conjugative Transposon" $tax_info ; then
				vir_name="Conjugative Transposon" ;
			elif grep -q "No homologues found" $tax_info ; then
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="genetic element" ;
				else
					vir_name="circular genetic element" ;
				fi
			elif grep -q "Circular genetic element" $tax_info ; then
				vir_name="Circular genetic element" ;
			elif grep -q "dsRNA virus" $tax_info ; then
				vir_name="dsRNA virus" ;
			elif grep -q "ssRNA virus" $tax_info ; then
				vir_name="ssRNA virus" ;
			elif grep -q "unclassified RNA virus" $tax_info ; then
				vir_name="unclassified RNA virus" ;
			elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
				vir_name="unclassified ssDNA bacterial virus" ;
			elif grep -q -i "phage" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "plasmid" $tax_info ; then
				vir_name="metagenomic plasmid" ;
			elif grep -q "Bacteria" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "unclassified virus" $tax_info ; then
				vir_name="virus" ;		
			elif grep -q -i "virus\|viridae" $tax_info ; then
				vir_name="virus" ;
			else
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="unclassified element" ;
				else
					vir_name="Circular genetic element" ;
				fi
			fi
		elif [ $TAX_CUTOFF == "order" ] ; then
			if grep -q "Virophage	Unclassified Taxon" $tax_info ; then
				vir_name="Virophage" ;
			elif grep -q "Adintovirus	Unclassified Taxon" $tax_info ; then
				vir_name="Adintovirus" ;
			elif grep -q "Polinton-like virus	Unclassified Taxon" $tax_info ; then
				vir_name="Polinton-like virus" ;
			elif grep -q "	order$" $tax_info ; then
				vir_name=$( grep "	order$" $tax_info | cut -f2 )
			elif grep -q "Conjugative Transposon" $tax_info ; then
				vir_name="Conjugative Transposon" ;
			elif grep -q "No homologues found" $tax_info ; then
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="genetic element" ;
				else
					vir_name="circular genetic element" ;
				fi
			elif grep -q "Circular genetic element" $tax_info ; then
				vir_name="Circular genetic element" ;
			elif grep -q "dsRNA virus" $tax_info ; then
				vir_name="dsRNA virus" ;
			elif grep -q "ssRNA virus" $tax_info ; then
				vir_name="ssRNA virus" ;
			elif grep -q "unclassified RNA virus" $tax_info ; then
				vir_name="unclassified RNA virus" ;
			elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
				vir_name="unclassified ssDNA bacterial virus" ;
			elif grep -q -i "phage" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "plasmid" $tax_info ; then
				vir_name="metagenomic plasmid" ;
			elif grep -q "Bacteria" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "unclassified virus" $tax_info ; then
				vir_name="virus" ;		
			elif grep -q -i "virus\|viridae" $tax_info ; then
				vir_name="virus" ;
			else
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="unclassified element" ;
				else
					vir_name="Circular genetic element" ;
				fi
			fi
		else
			if grep -q "Virophage	Unclassified Taxon" $tax_info ; then
				vir_name="Virophage" ;
			elif grep -q "Adintovirus	Unclassified Taxon" $tax_info ; then
				vir_name="Adintovirus" ;
			elif grep -q "Polinton-like virus	Unclassified Taxon" $tax_info ; then
				vir_name="Polinton-like virus" ;
			elif grep -q "Conjugative Transposon" $tax_info ; then
				vir_name="Conjugative Transposon" ;
			elif grep -q "No homologues found" $tax_info ; then
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="genetic element" ;
				else
					vir_name="circular genetic element" ;
				fi
			elif grep -q "Circular genetic element" $tax_info ; then
				vir_name="Circular genetic element" ;
			elif grep -q "dsRNA virus" $tax_info ; then
				vir_name="dsRNA virus" ;
			elif grep -q "ssRNA virus" $tax_info ; then
				vir_name="ssRNA virus" ;
			elif grep -q "unclassified RNA virus" $tax_info ; then
				vir_name="unclassified RNA virus" ;
			elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
				vir_name="unclassified ssDNA bacterial virus" ;
			elif grep -q -i "phage" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "plasmid" $tax_info ; then
				vir_name="metagenomic plasmid" ;
			elif grep -q "Bacteria" $tax_info ; then
				vir_name="Phage" ;
			elif grep -q "unclassified virus" $tax_info ; then
				vir_name="virus" ;		
			elif grep -q -i "virus\|viridae" $tax_info ; then
				vir_name="virus" ;
			else
				if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
					vir_name="unclassified element" ;
				else
					vir_name="Circular genetic element" ;
				fi
			fi
		fi

		#echo $vir_name ;
		fsa_head=$( echo $vir_name " sp." )
		tax_guess=$( tail -n+3 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out | grep "^[0-9]" | cut -f2 | tr '\n' ';' ) ; 
		perc_id=$( head -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out | sed 's/ /-/g' | awk '{FS="\t"; OFS="\t"} {print $2" "$3}' | sed 's/-/ /g' ) ;
		rand_id=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )

		# Editing and transferring tbl file and fasta (fsa) files to sequin directory

		if [ -s ${feat_tbl2%.comb3.tbl}.phan.fasta ]; then
			GCODE="11"
		else
			GCODE="1"
		fi
		if echo "$feat_tbl2" | grep -q "DTR_contigs_with_viral_domain" ; then
			
			NUCL_FILE="${feat_tbl2%.comb3.tbl}.rotate.fasta"
			input_contig_name=$( head -n1 ${feat_tbl2%.comb3.tbl}.rotate.fasta | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' )
			if [ "$WRAP" == "True" ] ; then 
				TOPOLOGY="circular"
			else
				TOPOLOGY="linear"
			fi
		else
			TOPOLOGY="linear"
			NUCL_FILE="${feat_tbl2%.comb3.tbl}.fna"
			if [ "$PROPHAGE" == "True" ] ; then
				if [ -s ${feat_tbl2%.comb3.tbl}.ITR.tbl ] ; then
					input_contig_name=$( head -n1 ${feat_tbl2%.comb3.tbl}.fna | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' )
				else
					input_contig_name=$( head -n1 ${feat_tbl2%_vs[0-9][0-9].comb3.tbl}.fna | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' )
				fi
			else
				input_contig_name=$( head -n1 ${feat_tbl2%.comb3.tbl}.fna | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' )
			fi
		fi
		INPUT_STEM=$( echo "$input_contig_name" | sed 's/@.*/@/g' )
		### crispr info
		if [ -s ${base_directory}/$CRISPR_FILE ] ; then
			if grep -q "$INPUT_STEM" ${base_directory}/$CRISPR_FILE ; then
				PRE_CRISP=$( grep -m1 "$INPUT_STEM" ${base_directory}/$CRISPR_FILE | awk '{OFS="\t"}{FS="\t"}{if ($3>=1) {print "yes"} else {print "no"}}' )
				if [ "$PRE_CRISP" == "yes" ] ;then
					CRISPR=$( grep -m1 "$INPUT_STEM" ${base_directory}/$CRISPR_FILE | cut -f2,3 | sed 's/	/: /g ; s/\(.*\)/\1 CRISPR spacer matches/'  )
				else 
					CRISPR=""
				fi
			else
				CRISPR=""
			fi
		else
			CRISPR=""
		fi
		if [ -s ${NUCL_FILE%.fna}.tax_guide.KNOWN_VIRUS.out ] ; then
			PRE_BLASTN=$( tail -n1 ${NUCL_FILE%.fna}.tax_guide.KNOWN_VIRUS.out )
			ANI=$( tail -n+2 ${NUCL_FILE%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${NUCL_FILE%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "; BLASTN hit: $PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
			fsa_head=$( echo $PRE_BLASTN " isolate" )
		elif [ -s ${NUCL_FILE%.fna}.tax_guide.CELLULAR.out ] ; then
			PRE_BLASTN=$( tail -n1 ${NUCL_FILE%.fna}.tax_guide.CELLULAR.out )
			ANI=$( tail -n+2 ${NUCL_FILE%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${NUCL_FILE%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "; BLASTN hit: $PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		elif [ -s ${NUCL_FILE%.rotate.fasta}.tax_guide.KNOWN_VIRUS.out ] ; then
			PRE_BLASTN=$( tail -n1 ${NUCL_FILE%.rotate.fasta}.tax_guide.KNOWN_VIRUS.out )
			ANI=$( tail -n+2 ${NUCL_FILE%.rotate.fasta}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${NUCL_FILE%.rotate.fasta}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "; BLASTN hit: $PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
			fsa_head=$( echo $PRE_BLASTN " isolate" )
		elif [ -s ${NUCL_FILE%.rotate.fasta}.tax_guide.CELLULAR.out ] ; then
			PRE_BLASTN=$( tail -n1 ${NUCL_FILE%.rotate.fasta}.tax_guide.CELLULAR.out )
			ANI=$( tail -n+2 ${NUCL_FILE%.rotate.fasta}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${NUCL_FILE%.rotate.fasta}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "; BLASTN hit: $PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		else
			BLASTN_INFO=""
		fi

		cp $feat_tbl2 sequin_and_genome_maps/${JUST_TBL2_FILE%.comb3.tbl}.tbl ; 
				#echo "1"
		bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -v topoq="$TOPOLOGY" -v gcodeq="$GCODE" -v o_name="$input_contig_name" -v crispr1="$CRISPR" -v blastn="$BLASTN_INFO" -c fastx '{ print ">" newname " [note=input name:"o_name" -- closest relative: " tax_var " " perc_var " ; " crispr1" "blastn"] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology="topoq"] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode="gcodeq"]" ; print $seq }' $NUCL_FILE > sequin_and_genome_maps/${JUST_TBL2_FILE%.comb3.tbl}.fsa ; 

	#echo $input_contig_name
	if [ -s reads_to_all_contigs_over${LENGTH_MINIMUM}nt.coverage.txt ] ; then
		if echo "$JUST_TBL2_FILE" | grep -q "_vs[0-9][0-9]" ; then
			COVERAGE=$( grep "${JUST_TBL2_FILE%_vs[0-9][0-9].comb3.tbl}	" reads_to_all_contigs_over${LENGTH_MINIMUM}nt.coverage.txt | cut -f2 )
		else
			COVERAGE=$( grep "${JUST_TBL2_FILE%.comb3.tbl}	" reads_to_all_contigs_over${LENGTH_MINIMUM}nt.coverage.txt | cut -f2 )
		fi

		#echo $COVERAGE
	else
		COVERAGE="1"
	fi
	echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > sequin_and_genome_maps/${file_core}.cmt ;
	echo "Assembly Method	" $ASSEMBLER >> sequin_and_genome_maps/${file_core}.cmt ;
	echo "Genome Coverage	"$COVERAGE"x" >> sequin_and_genome_maps/${file_core}.cmt ;
	echo "Sequencing Technology	Illumina" >> sequin_and_genome_maps/${file_core}.cmt ;
	echo "Annotation Pipeline	Cenote-Taker2" >> sequin_and_genome_maps/${file_core}.cmt ;
	echo "URL	https://github.com/mtisza1/Cenote-Taker2" >> sequin_and_genome_maps/${file_core}.cmt ;		
	done
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running tbl2asn " $MDYT
	if [[ $DATA_SOURCE = "tpa_assembly" ]] ;then
		tbl2asn -V vb -j "[keyword=TPA:assembly]" -t ${base_directory}/${template_file} -X C -p sequin_and_genome_maps/ ;
	else
		tbl2asn -V vb -t ${base_directory}/${template_file} -X C -p sequin_and_genome_maps/ ;
	fi
else
	echo "no tbl file found for sequin/genome map"
fi

# make gtf tables
echo "Making gtf tables from final feature tables"
if [ -d sequin_and_genome_maps ] ; then
	cd sequin_and_genome_maps
	COMB4_TBL=$( find . -maxdepth 1 -type f -name "*.tbl" )
	if [ -n "$COMB4_TBL" ] ; then
		for feat_tbl2 in $COMB4_TBL ; do
			if [ -s ${feat_tbl2%.tbl}.gtf ] ; then
				rm -f ${feat_tbl2%.tbl}.gtf
			fi
			grep "^[0-9]\|^<[0-9]" -A3 $feat_tbl2 | sed '/--/d' | sed 's/ /_/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n		//g' | while read LINE ; do
				if echo $LINE | grep -q "CDS" ; then
					GENOME=${feat_tbl2%.tbl}
					FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
					FEAT_START=$( echo $LINE | cut -d " " -f1 )
					FEAT_END=$( echo $LINE | cut -d " " -f2 )
					FEAT_NAME=$( echo $LINE | cut -d " " -f7 )
					FEAT_ATT=$( echo $LINE | cut -d " " -f9 )
					FEAT_ID=$( echo $LINE | cut -d " " -f5 )
				elif echo $LINE | grep -q "repeat_region" && grep -q "DTR" ; then
					GENOME=${feat_tbl2%.tbl}
					FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
					FEAT_START=$( echo $LINE | cut -d " " -f1 )
					FEAT_END=$( echo $LINE | cut -d " " -f2 )
					FEAT_NAME="DTR"
					FEAT_ATT="DTR"
					FEAT_ID="DTR"	
				elif echo $LINE | grep -q "repeat_region" ; then
					GENOME=${feat_tbl2%.tbl}
					FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
					FEAT_START=$( echo $LINE | cut -d " " -f1 )
					FEAT_END=$( echo $LINE | cut -d " " -f2 )
					FEAT_NAME="ITR"
					FEAT_ATT="ITR"
					FEAT_ID="ITR"	
				fi
				echo -e "$GENOME\t""Cenote-Taker\t""$FEAT_TYPE\t""$FEAT_START\t""$FEAT_END\t"".\t"".\t"".\t""gene_id \"$FEAT_ID\"; gene_name \"$FEAT_NAME\"; gene_inference \"$FEAT_ATT\"" >> ${feat_tbl2%.tbl}.gtf
			done
		done
	fi
	cd ${base_directory}/${run_title}
fi
### conjugative machinery table
if [ -d sequin_and_genome_maps ] ; then
	cd sequin_and_genome_maps
	COMB4_TBL=$( find . -maxdepth 1 -type f -name "*.tbl" )
	if [ -n "$COMB4_TBL" ] ; then
		for feat_tbl2 in $COMB4_TBL ; do
			if [ -s ${feat_tbl2%.tbl}.gtf ] ; then
				CONJ_COUNT=$( grep -i "virb\|type-IV\|secretion system\|conjugation\|conjugal\|transposon\| tra[a-z] \|	tra[a-z]\|trb[b-z]\|pilus" $feat_tbl2 | grep -v "TRAF\|TRAP\|protein tyrosine phosphatase\|ttRBP\|SpoU\|transport\|central" | wc -l )
				STRUCTURAL_COUNT=$( grep -i "capsid\|terminase\|portal\|baseplate\|base plate\|tail\|collar\|zot\|zonular\|minor coat\|packaging\|	virion protein" $feat_tbl2 | wc -l )
				if [[ $CONJ_COUNT -gt 0 ]] && [[ $STRUCTURAL_COUNT == 0 ]] ; then
					grep -v "TRAF\|TRAP\|protein tyrosine phosphatase\|ttRBP\|SpoU\|transport" $feat_tbl2 | grep -B2 -i "virb\|type-IV\|secretion system\|conjugation\|conjugal\|transposon\| tra[a-z] \|	tra[a-z]\|trb[b-z]\|pilus" | grep "^[0-9]\|^<[0-9]" | cut -f1,2 | while read START_STOP ; do 
						grep "$START_STOP" ${feat_tbl2%.tbl}.gtf
					done > ${feat_tbl2%.tbl}.putative_conjugative_machinery.gtf
				fi
			fi
		done
	fi
	cd ${base_directory}/${run_title}
fi
###

echo " "
MDYT=$( date +"%m-%d-%y---%T" )
echo "time update: Finishing " $MDYT
if [ -s final_combined_virus_sequences_${run_title}.fna ] ; then
	NUMBER_VIRUSES=$( grep -F ">" final_combined_virus_sequences_${run_title}.fna | wc -l | bc )
	DTR_TOTAL=$( grep -F " DTR" final_combined_virus_sequences_${run_title}.fna | wc -l | bc )
	ITR_TOTAL=$( grep -F " ITR" final_combined_virus_sequences_${run_title}.fna | wc -l | bc )
	NO_END_TOTAL=$( grep -F " no_end_feature" final_combined_virus_sequences_${run_title}.fna | wc -l | bc )
	echo "$(tput setaf 3)Virus prediction summary:$(tput sgr 0)"
	echo "$NUMBER_VIRUSES virus contigs were detected/predicted. $DTR_TOTAL contigs had DTRs/circularity. $ITR_TOTAL contigs had ITRs. $NO_END_TOTAL were linear/had no end features"
fi
if [ -s ${run_title}_PRUNING_INFO_TABLE.tsv ] ; then
	PRUNE_ATTEMPTS=$( awk '{OFS="\t"}{FS="\t"}{ if ($6=="True") {print}}' ${run_title}_PRUNING_INFO_TABLE.tsv | wc -l | bc )
	PRUNE_REMOVE=$( awk '{OFS="\t"}{FS="\t"}{ if ($6=="True" && $7=="True") {print}}' ${run_title}_PRUNING_INFO_TABLE.tsv | wc -l | bc )
	PRUNE_REMAIN=$( awk '{OFS="\t"}{FS="\t"}{ if ($6=="True" && $7=="False") {print}}' ${run_title}_PRUNING_INFO_TABLE.tsv | wc -l | bc )
	echo "$(tput setaf 3)Prophage pruning summary:$(tput sgr 0)"
	echo "$PRUNE_ATTEMPTS linear contigs > 10 kb were run through pruning module, and $PRUNE_REMOVE virus sub-contigs (putative prophages/proviruses) were extracted from these. $PRUNE_REMAIN virus contigs were kept intact."

fi

echo "ORIGINAL_NAME	CENOTE_NAME	ORGANISM_NAME	END_FEATURE	LENGTH	ORF_CALLER	NUM_HALLMARKS	HALLMARK_NAMES	BLASTN_INFO" > ${run_title}_CONTIG_SUMMARY.tsv
CIRCULAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -wholename "DTR_contigs_with_viral_domain/*fna" | grep -v "\.rc\.fna" )

if [ -n "$CIRCULAR_HALLMARK_CONTIGS" ] ; then
	for LINEAR in $CIRCULAR_HALLMARK_CONTIGS ; do 
		CENOTE_NAME=$( head -n1 $LINEAR | cut -d " " -f1 | sed 's/>//g' )
		ORIGINAL_NAME=$( head -n1 $LINEAR | cut -d " " -f2 )
		LENGTH=$( bioawk -c fastx '{print length($seq)}' $LINEAR )
		if [ -s ${LINEAR%.fna}.rotate.AA.hmmscan.sort.out ] ; then
			NUM_HALLMARKS=$( cat ${LINEAR%.fna}.rotate.AA.hmmscan.sort.out | wc -l | bc )
			HALLMARK_NAMES=$( cut -f1 ${LINEAR%.fna}.rotate.AA.hmmscan.sort.out | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' | sort -u | sed 's/,//g' )
		else
			NUM_HALLMARKS=0
			HALLMARK_NAMES="none"
		fi
		JUST_FILE=$( echo "$LINEAR" | sed 's/.*\///g' )
		if [ -s DTR_contigs_with_viral_domain/DTR_seqs_for_phanotate.txt ] ; then
			if grep -q "${JUST_FILE%.fna}.rotate.fasta" DTR_contigs_with_viral_domain/DTR_seqs_for_phanotate.txt ; then
				ORF_CALLER="PHANOTATE"
			else
				ORF_CALLER="Prodigal"
			fi
		else
			ORF_CALLER="Prodigal"
		fi
		END_FEATURE="DTR"
		if [ -s sequin_and_genome_maps/${CENOTE_NAME}.fsa ] ; then
			ORGANISM=$( head -n1 sequin_and_genome_maps/${CENOTE_NAME}.fsa | sed 's/.*\[organism=\(.*\)\] \[moltype=.*/\1/' )
		else
			ORGANISM="Unknown"
		fi
		if [ -s ${LINEAR%.fna}.tax_guide.KNOWN_VIRUS.out ] ; then
			PRE_BLASTN=$( tail -n1 ${LINEAR%.fna}.tax_guide.KNOWN_VIRUS.out )
			ANI=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "$PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		elif [ -s ${LINEAR%.fna}.tax_guide.CELLULAR.out ] ; then
			PRE_BLASTN=$( tail -n1 ${LINEAR%.fna}.tax_guide.CELLULAR.out )
			ANI=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "$PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		elif [ "$handle_knowns" != "blast_knowns" ] ; then
			BLASTN_INFO="BLASTN not conducted"
		else
			BLASTN_INFO="no high coverage hits"
		fi		
		echo "${ORIGINAL_NAME}	${CENOTE_NAME}	${ORGANISM}	${END_FEATURE}	${LENGTH}	${ORF_CALLER}	${NUM_HALLMARKS}	${HALLMARK_NAMES}	${BLASTN_INFO}" >> ${run_title}_CONTIG_SUMMARY.tsv
	done
fi

if [ "$PROPHAGE" == "True" ] ; then
	LIST_OF_VIRAL_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "no_end_contigs_with_viral_domain/*_vs[0-9][0-9].fna" )
else
	LIST_OF_VIRAL_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "no_end_contigs_with_viral_domain/*fna" )
fi

if [ -n "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ; then
	for LINEAR in $LIST_OF_VIRAL_DOMAIN_CONTIGS ; do 
		CENOTE_NAME=$( head -n1 $LINEAR | cut -d " " -f1 | sed 's/>//g' )
		if [ "$PROPHAGE" == "True" ] ; then
			if [ -s ${LINEAR%.fna}.ITR.tbl ] ; then
				ORIGINAL_NAME=$( head -n1 $LINEAR | cut -d " " -f2 )
			else
				ORIGINAL_NAME=$( head -n1 ${LINEAR%_vs[0-9][0-9].fna}.fna | cut -d " " -f2 )
			fi
		else
			ORIGINAL_NAME=$( head -n1 $LINEAR | cut -d " " -f2 )
		fi

		LENGTH=$( bioawk -c fastx '{print length($seq)}' $LINEAR )
		if [ -s ${LINEAR%.fna}.AA.hmmscan.sort.out ] ; then
			NUM_HALLMARKS=$( cat ${LINEAR%.fna}.AA.hmmscan.sort.out | wc -l | bc )
			HALLMARK_NAMES=$( cut -f1 ${LINEAR%.fna}.AA.hmmscan.sort.out | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' | sort -u | sed 's/,//g' )
		else
			NUM_HALLMARKS=0
			HALLMARK_NAMES="none"
		fi
		if [ -s ${LINEAR%.fna}.ITR.tbl ] ; then
			END_FEATURE="ITR"
		else
			END_FEATURE="None"
		fi
		JUST_FILE=$( echo "$LINEAR" | sed 's/.*\///g' )
		if [ -s no_end_contigs_with_viral_domain/LIN_seqs_for_phanotate.txt ] ; then 
			if grep -q "${JUST_FILE}" no_end_contigs_with_viral_domain/LIN_seqs_for_phanotate.txt ; then
				ORF_CALLER="PHANOTATE"
			else
				ORF_CALLER="Prodigal"
			fi
		else
			ORF_CALLER="Prodigal"
		fi
		if [ -s sequin_and_genome_maps/${CENOTE_NAME}.fsa ] ; then
			ORGANISM=$( head -n1 sequin_and_genome_maps/${CENOTE_NAME}.fsa | sed 's/.*\[organism=\(.*\)\] \[moltype=.*/\1/' )
		else
			ORGANISM="Unknown"
		fi
		if [ -s ${LINEAR%.fna}.tax_guide.KNOWN_VIRUS.out ] ; then
			PRE_BLASTN=$( tail -n1 ${LINEAR%.fna}.tax_guide.KNOWN_VIRUS.out )
			ANI=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "$PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		elif [ -s ${LINEAR%.fna}.tax_guide.CELLULAR.out ] ; then
			PRE_BLASTN=$( tail -n1 ${LINEAR%.fna}.tax_guide.CELLULAR.out )
			ANI=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f3 )
			AF=$( tail -n+2 ${LINEAR%.fna}.blastn_intraspecific.out | head -n1 | cut -f4 )
			BLASTN_INFO=$( echo "$PRE_BLASTN", ANI=${ANI}%, AF=${AF}% )
		elif [ "$handle_knowns" != "blast_knowns" ] ; then
			BLASTN_INFO="BLASTN not conducted"
		else
			BLASTN_INFO="no high coverage hits"
		fi	
		echo "${ORIGINAL_NAME}	${CENOTE_NAME}	${ORGANISM}	${END_FEATURE}	${LENGTH}	${ORF_CALLER}	${NUM_HALLMARKS}	${HALLMARK_NAMES}	${BLASTN_INFO}" >> ${run_title}_CONTIG_SUMMARY.tsv
	done
fi

echo "removing ancillary files"

if [ -d DTR_contigs_with_viral_domain ] ; then
	cd DTR_contigs_with_viral_domain
	rm -f *.all_start_stop.txt *.bad_starts.txt *.comb.tbl *.comb2.tbl *.good_start_orfs.txt *.hypo_start_stop.txt *.nucl_orfs.fa *.remove_hypo.txt *.log *.promer.contigs_with_ends.fa *.promer.promer *.out.hhr *.starting_orf.txt *.out.hhr *.nucl_orfs.txt *.called_hmmscan.txt *.hmmscan_replicate.out *.hmmscan.out *.rotate.no_hmmscan.fasta *.starting_orf.1.fa *.phan.*fasta *used_positions.txt *.prodigal.for_prodigal.fa *.prodigal.gff *.trnascan-se2.txt *.for_blastp.txt *.for_hhpred.txt circular_contigs_spades_names.txt SPLIT_CIRCULAR_AA*fasta all_circular_contigs_${run_title}.fna SPLIT_DTR_* *called_hmmscan*txt *HH.tbl *CDS.bed *tRNA.bed *ORFs_over_tRNAs.tsv *prodigal.fasta *blast_hypo.fasta *no_hmmscan2.fasta *no_hmmscan1.fasta *all_hhpred_queries.AA.fasta
	cd ${base_directory}/${run_title}
fi
rm -rf bt2_indices/

find other_contigs/ -type f -name "*.AA.fasta" -exec rm -f {} \;
find other_contigs/ -type f -name "*.AA.sorted.fasta" -exec rm -f {} \;
find other_contigs/ -type f -name "*.out" -exec rm -f {} \;
find other_contigs/ -type f -name "*.dat" -exec rm -f {} \;
find other_contigs/ -type f -name "*.called_hmmscan.txt" -exec rm -f {} \;
find other_contigs/ -type f -name "SPLIT_LARGE_GENOME_AA_*fasta" -exec rm -f {} \;
find ITR_containing_contigs/ -type f -name "SPLIT_ITR_AA*fasta" -exec rm -f {} \;
find . -type f -name "SPLIT_CIRCULAR_AA*" -exec rm -f {} \;
rm -f *called_hmmscan.txt circular_contigs_spades_names.txt
if [ -d no_end_contigs_with_viral_domain ] ; then
	find no_end_contigs_with_viral_domain/ -type f -name "*.called_hmmscan2.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*.hmmscan2.out" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*all_hhpred_queries.AA.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*all_start_stop.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*trnascan-se2.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*for_hhpred.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*for_blastp.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*HH.tbl" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*hypo_start_stop.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*remove_hypo.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*rps_nohits.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*tax_guide.blastx.tab" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*trans.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*called_hmmscan*.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*no_hmmscan*.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "SPLIT_LIN_HMM2_GENOME_AA*fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "SPLIT_LIN_sort_GENOME_AA*" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "SPLIT_LIN_RPS_AA*" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*used_positions.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*seq_chunk_coordinates.csv" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*blast_hypo.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*.bed" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*ORFs_over_tRNAs.tsv" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*prodigal.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*all_called_hmmscans.txt" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*phan.fasta" -exec rm -f {} \;
	find no_end_contigs_with_viral_domain/ -type f -name "*phan.sort.fasta" -exec rm -f {} \;
fi

echo "$(tput setaf 3)output directory: "$run_title" $(tput sgr 0)"
echo "$(tput setaf 3) >>>>>>CENOTE-TAKER 2 HAS FINISHED TAKING CENOTES<<<<<< $(tput sgr 0)"



