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
echo "Version 2.1.2"
echo " "

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

if [ "$ANNOTATION_MODE" == "True" ] ; then
	LIN_MINIMUM_DOMAINS=0
	CIRC_MINIMUM_DOMAINS=0
	circ_length_cutoff=1000
	linear_length_cutoff=1
fi

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
echo "Annotation mode?                   $ANNOTATION_MODE"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@"

if [ "${SCRATCH_DIR}" == "none" ] ; then
	echo "scratch space will not be used in this run"
	CD_HHSUITE="${CENOTE_SCRIPT_DIR}/NCBI_CD/NCBI_CD"
	PFAM_HHSUITE="${CENOTE_SCRIPT_DIR}/pfam_32_db/pfam"
	PDB_HHSUITE="${CENOTE_SCRIPT_DIR}/pdb70/pdb70"
	echo "HHsuite database locations:"
	echo $CD_HHSUITE
	echo $PFAM_HHSUITE
	echo $PDB_HHSUITE
elif [ -d ${SCRATCH_DIR}/ ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: setting up lscratch databases: " $MDYT
	if [ ! -s ${SCRATCH_DIR}/NCBI_CD/NCBI_CD_a3m.ffdata ] ; then
		mkdir ${SCRATCH_DIR}/NCBI_CD
		cp ${CENOTE_SCRIPT_DIR}/NCBI_CD/NCBI_CD* ${SCRATCH_DIR}/NCBI_CD/
	fi
	if [ ! -s ${SCRATCH_DIR}/pfam_32_db/pfam_a3m.ffdata ] ; then	
		mkdir ${SCRATCH_DIR}/pfam_32_db
		cp ${CENOTE_SCRIPT_DIR}/pfam_32_db/pfam* ${SCRATCH_DIR}/pfam_32_db/
	fi
	if [ ! -s ${SCRATCH_DIR}/pdb70/pdb70_a3m.ffdata ] ; then		
		mkdir ${SCRATCH_DIR}/pdb70
		cp ${CENOTE_SCRIPT_DIR}/pdb70/pdb70* ${SCRATCH_DIR}/pdb70/
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

if [[ ":$PATH:" == *":${CENOTE_SCRIPT_DIR}/hh-suite/build/scripts:"* ]] && [[ ":$PATH:" == *":${CENOTE_SCRIPT_DIR}/hh-suite/build/bin:"* ]] ; then 
	echo "hhsuite is in path" 
else 
	export PATH="${CENOTE_SCRIPT_DIR}/hh-suite/build/bin:${CENOTE_SCRIPT_DIR}/hh-suite/build/scripts:$PATH" 
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


if [ ${original_contigs: -6} == ".fasta" ]; then
	echo "$(tput setaf 5)File with .fasta extension detected, attempting to keep contigs over $LENGTH_MINIMUM nt and find circular sequences with apc.pl$(tput sgr 0)"
	bioawk -v run_var="$run_title" -v contig_cutoff="$LENGTH_MINIMUM" -c fastx '{ if(length($seq) > contig_cutoff) { print ">"run_var NR" "$name; print $seq }}' $original_contigs > ${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta ;
	cd $run_title
	echo "cenote_shortcut" > ${run_title}_CONTIG_SUMMARY.tsv
	perl ${CENOTE_SCRIPT_DIR}/apc_cenote1.pl -b $run_title -c $CENOTE_SCRIPT_DIR ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta >/dev/null 2>&1
	rm -f apc_aln*
	APC_CIRCS=$( find * -maxdepth 0 -type f -name "${run_title}*.fa" )
	if [ -n "$APC_CIRCS" ] ;then
		for fa1 in $APC_CIRCS ; do
			CIRC_SEQ_NAME=$( head -n1 $fa1 | sed 's/|.*//g' ) ; 
			CIRC_NEW_NAME=$( echo "$CIRC_SEQ_NAME" | sed 's/>//g ; s/ .*//g' )
			grep -A1 "^$CIRC_SEQ_NAME" ../${original_contigs%.fasta}.over_${LENGTH_MINIMUM}nt.fasta | sed '/--/d' > ${CIRC_NEW_NAME}.fasta
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
	perl ${CENOTE_SCRIPT_DIR}/apc_cenote1.pl -b $run_title -c $CENOTE_SCRIPT_DIR ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta >/dev/null 2>&1
	rm -f apc_aln*
	APC_CIRCS=$( find * -maxdepth 0 * -type f -name "${run_title}*.fa" )
	if [ -n "$APC_CIRCS" ] ;then
		for fa1 in $APC_CIRCS ; do 
			CIRC_SEQ_NAME=$( head -n1 $fa1 | sed 's/|.*//g' ) ; 
			CIRC_NEW_NAME=$( echo "$CIRC_SEQ_NAME" | sed 's/>//g ; s/ .*//g' )
			grep -A1 "^$CIRC_SEQ_NAME" ../${original_contigs%.fastg}.over_${LENGTH_MINIMUM}nt.fasta | sed '/--/d' > ${CIRC_NEW_NAME}.fasta
			echo "${CIRC_NEW_NAME}.fasta is circular"
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
CIRC_CONTIGS=$( find * -maxdepth 0 -type f -name "*.fasta" )
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
original_fastas=$( find * -maxdepth 0 -type f -name "*.fasta" )
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
CONTIGS_NON_CIRCULAR=$( find * -maxdepth 0 -type f -name "*[0-9].fasta" )

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
						echo ${DAT%.2.3.5.80.10.40.500000.10000.dat} "contains ITRs:" ; 
						echo $LENGTH "5-prime ITR:" $LOW_START "3-prime ITR:" $HIGH_END ; 
						mv ${DAT%.2.3.5.80.10.40.500000.10000.dat} ../ITR_containing_contigs/${DAT%.fasta.2.3.5.80.10.40.500000.10000.dat}.fasta
						echo "$(tput setaf 4) Making ITR .tbl file $(tput sgr 0)"
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
	CONTIGS_NON_CIRCULAR=$( find * -maxdepth 0 -type f -name "*[0-9].fasta" )
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

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running linear contigs with hmmscan against virus hallmark gene database: $virus_domain_db " $MDYT	
		cat $( find * -maxdepth 0 -type f -name "*.AA.sorted.fasta" ) > all_large_genome_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_large_genome_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_LARGE_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_large_genome_proteins.AA.fasta
			SPLIT_AA_LARGE=$( find * -maxdepth 0 -type f -name "SPLIT_LARGE_GENOME_AA_*.fasta" )
			if  [[ $virus_domain_db = "standard" ]] ; then
				echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
				echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
			elif [[ $virus_domain_db = "virion" ]]; then
				echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
			elif [[ $virus_domain_db = "rna_virus" ]]; then
				echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
			else
				echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
				rm -f ./*{0..9}.fasta
				break
			fi		
			HMM_REP_NUMEBR=$( find * -maxdepth 0 -type f -name "SPLIT_LARGE_GENOME_AA_*AA.hmmscan_replicate.out" | wc -l )
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
			if [ -s LARGE_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
				cut -f3 LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
					HALL_COUNT=$( grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | wc -l | bc )
					if [ $HALL_COUNT -ge $LIN_MINIMUM_DOMAINS ] ; then 
						mv ${HIT}.fasta ../no_end_contigs_with_viral_domain/${HIT}.fna
						grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out > ../no_end_contigs_with_viral_domain/${HIT}.AA.hmmscan.sort.out
						# ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.no_hmmscan1.fasta
						grep "${HIT}_" LARGE_GENOME_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan.txt
						grep -v -f ${HIT}.AA.called_hmmscan.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ../no_end_contigs_with_viral_domain/${HIT}.AA.no_hmmscan1.fasta
						mv ${HIT}.AA.sorted.fasta ../no_end_contigs_with_viral_domain/

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
				cat $REMAINDER >> non_viral_domains_contigs.fna
				rm -f $REMAINDER
			fi
		done
	fi
fi


cd ..


DTR_SEQS=$( find * -maxdepth 0 -type f -regextype sed -regex "${run_title}[0-9]\{1,6\}.fasta" )


if [ ! -z "$DTR_SEQS" ] ; then
	mkdir DTR_contigs_with_viral_domain
	if [ $CIRC_MINIMUM_DOMAINS -le 0 ] ; then
		for CIRC in $DTR_SEQS ; do
			mv ${CIRC} DTR_contigs_with_viral_domain/${CIRC%.fasta}.fna
		done
	else
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: Calling ORFs for circular/DTR sequences with prodigal " $MDYT
		echo "$DTR_SEQS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t prodigal -a {}.AA.fasta -i {}.fasta -p meta -q >/dev/null 2>&1
		for CIRC in $DTR_SEQS ; do 
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
		cat $( find * -maxdepth 0 -type f -name "*.AA.sorted.fasta" ) > all_circular_genome_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_circular_genome_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_CIRCULAR_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_circular_genome_proteins.AA.fasta
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running hmmscan on circular/DTR contigs " $MDYT
			SPLIT_AA_CIRC=$( find * -maxdepth 0 -type f -name "SPLIT_CIRCULAR_AA_*.fasta" )
			if  [[ $virus_domain_db = "standard" ]] ; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_replication_clusters3 {}.fasta	>/dev/null 2>&1
			elif [[ $virus_domain_db = "virion" ]]; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta	>/dev/null 2>&1	
			elif [[ $virus_domain_db = "rna_virus" ]]; then
				echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
			else
				echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
				rm -f ./*{0..9}.fasta
				break
			fi
			HMM_REP_NUMEBR=$( find * -maxdepth 0 -type f -name "SPLIT_CIRCULAR_AA_*AA.hmmscan_replicate.out" | wc -l )
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
						rm ${HIT}.AA.fasta

					elif [ -s ${HIT}.fasta ] ; then
						sed 's/ /#/g' ${HIT}.fasta | bioawk -c fastx '{print ">"$name"#DTRs" ; print $seq}' | 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
						rm -f ${HIT}.fasta
					fi
				done
			fi
		else
			echo "No AA seqs found in circular contigs, discover viruses module"
		fi
		for REMAINDER in $DTR_SEQS ; do
			if [ -s $REMAINDER ] ; then
				sed 's/ /#/g' $REMAINDER | bioawk -c fastx '{print ">"$name"#DTRs" ; print $seq}' | 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
				rm -f $REMAINDER
			fi
		done
	fi
fi

ITR_SEQS=$( find * -maxdepth 1 -type f -wholename "ITR_containing_contigs/*fasta" )


if [ ! -z "$ITR_SEQS" ] ; then 
	cd ITR_containing_contigs
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
	cat $( find * -maxdepth 0 -type f -name "*.AA.sorted.fasta" ) > all_ITR_genome_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_ITR_genome_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_ITR_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_ITR_genome_proteins.AA.fasta
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running hmmscan on ITR contigs " $MDYT
		SPLIT_AA_CIRC=$( find * -maxdepth 0 -type f -name "SPLIT_ITR_AA_*.fasta" )
		if  [[ $virus_domain_db = "standard" ]] ; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
		elif [[ $virus_domain_db = "virion" ]]; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
		elif [[ $virus_domain_db = "rna_virus" ]]; then
			echo "$SPLIT_AA_CIRC" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
		else
			echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
			rm -f ./*{0..9}.fasta
			break
		fi
		HMM_REP_NUMEBR=$( find * -maxdepth 0 -type f -name "SPLIT_ITR_AA_*AA.hmmscan_replicate.out" | wc -l )
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
		if [ -s ITR_GENOME_COMBINED.AA.hmmscan.sort.out ] ; then
			cut -f3 ITR_GENOME_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
				HALL_COUNT=$( grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out | wc -l | bc )
				if [ $HALL_COUNT -ge $CIRC_MINIMUM_DOMAINS ] ; then 
					mv ${HIT}.fasta ${HIT}.fna
					grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out > ${HIT}.AA.hmmscan.sort.out
					grep "${HIT}_" ITR_GENOME_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan.txt
					grep -v -f ${HIT}.AA.called_hmmscan.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.AA.no_hmmscan1.fasta
					rm ${HIT}.AA.fasta

				elif [ -s ${HIT}.fasta ] ; then
					sed 's/ /#/g' ${HIT}.fasta | bioawk -c fastx '{print ">"$name"#ITRs" ; print $seq}' | 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
					rm -f ${HIT}.fasta
				fi
			done
		fi
	else
		echo "no AA seqs found in ITR contigs, discover viruses module"
	fi

	for REMAINDER in $ITR_SEQS ; do
		if [ -s $REMAINDER ] ; then
			sed 's/ /#/g' $REMAINDER | bioawk -c fastx '{print ">"$name"#ITRs" ; print $seq}' | 's/#/ /g' >> other_contigs/non_viral_domains_contigs.fna
			rm -f $REMAINDER
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
CIRCULAR_HALLMARK_CONTIGS=$( find * -maxdepth 0 -type f -name "*fna" )
if [ -n "$CIRCULAR_HALLMARK_CONTIGS" ] ; then
	echo "Annotating DTR contigs"
	#echo $PWD
	for nucl_fa in $CIRCULAR_HALLMARK_CONTIGS ; do
		#echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"
		getorf -circular -minsize 240 -table 11 -find 3 -sequence $nucl_fa -outseq ${nucl_fa%.fna}.nucl_orfs.fa >/dev/null 2>&1 

		grep ">" ${nucl_fa%.fna}.nucl_orfs.fa > ${nucl_fa%.fna}.nucl_orfs.txt
		cat "${nucl_fa%.fna}.nucl_orfs.txt" | while read liner ; do
			start_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
			end_base=$( echo $liner | sed 's/.*- \(.*\)\].*/\1/' )
			length=$(( $start_base-$end_base ))
			abso_length=$( echo $length | sed 's/-//g' )
			if [ $abso_length -gt 239 ]; then
				if [[ "$end_base" -gt "$start_base" ]]; then
					for ((counter_f=(( $start_base + 1 ));counter_f<=(( $end_base + 3 ));counter_f++)); do
						echo "$counter_f" >> ${nucl_fa%.fna}.bad_starts.txt
						
					done
				elif [[ "$start_base" -gt "$end_base" ]]; then
					for ((counter_r=(( $end_base - 3 ));counter_r<=(( $start_base - 1 ));counter_r++)) ; do
						echo "$counter_r" >> ${nucl_fa%.fna}.bad_starts.txt
					done
				fi
			fi
		done

		cat "${nucl_fa%.fna}.nucl_orfs.txt" | while read liner ; do
			starter_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
			if grep -q "$starter_base" ${nucl_fa%.fna}.bad_starts.txt ; then
				continue
			else
				echo $liner >> ${nucl_fa%.fna}.good_start_orfs.txt
			fi
		done
		if [ -s "${nucl_fa%.fna}.good_start_orfs.txt" ]; then	
			cut -d " " -f1 ${nucl_fa%.fna}.good_start_orfs.txt | head -n1 | sed 's/>//g' > ${nucl_fa%.fna}.starting_orf.txt
			bioawk -c fastx '{ print $name, $seq, length($seq) }' ${nucl_fa%.fna}.nucl_orfs.fa | grep -f ${nucl_fa%.fna}.starting_orf.txt | sed '/--/d' | head -n1 | awk '{print ">"$1, $3; print $2}' > ${nucl_fa%.fna}.starting_orf.1.fa ;
			circlator fixstart --genes_fa ${nucl_fa%.fna}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fna}.rotate ;
		else
			head -n2 ${nucl_fa%.fna}.nucl_orfs.fa | sed '/--/d' > ${nucl_fa%.fna}.starting_orf.1.fa ;
			circlator fixstart --genes_fa ${nucl_fa%.fna}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fna}.rotate ;
		fi
	done
fi

#-# blastx for translation decision
ROTATED_DTR_CONTIGS=$( find * -maxdepth 0 -type f -regextype sed -regex "${run_title}[0-9]\{1,6\}.rotate.fasta" )

if [ -n "$ROTATED_DTR_CONTIGS" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running BLASTX, DTR contigs " $MDYT
	echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU -t blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -num_threads 1 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query {}.rotate.fasta -out {}.tax_guide.blastx.out >/dev/null 2>&1
	echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta/.fasta/g' | while read nucl_fa ; do
		if [ ! -s "${nucl_fa%.fasta}.tax_guide.blastx.out" ]; then
			echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out ;
		elif grep -i -q "circovir\|genomovir\|geminivir\|nanovir\|redondovir\|bacilladnavir\|smacovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then 
			EVALUE=$( head -n1 "${nucl_fa%.fasta}.tax_guide.blastx.out" | cut -f4 ) ; 
			NEW_TAX=$( head -n1 ${nucl_fa%.fasta}.tax_guide.blastx.out | awk -v VALUE="$EVALUE" '{if (VALUE>1e-50) { print $0 ; print "CRESS virus" } else { print $0}}' )
			echo "$NEW_TAX" > ${nucl_fa%.fasta}.tax_guide.blastx.out ;
			if grep -q "CRESS virus" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
				echo ${nucl_fa%.fasta} "is a CRESS virus"
			else
				ktClassifyBLAST -o ${nucl_fa%.fasta}.tax_guide.blastx.tab ${nucl_fa%.fasta}.tax_guide.blastx.out >/dev/null 2>&1
				taxid=$( tail -n1 ${nucl_fa%.fasta}.tax_guide.blastx.tab | cut -f2 )
				efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element Lineage >> ${nucl_fa%.fasta}.tax_guide.blastx.out
				sleep 0.4s
			fi
		elif grep -q "virophage" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Virophage" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		elif grep -q "adinto" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Adintovirus" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		elif grep -i -q "polinton" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo "Polinton-like virus" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
		else
			ktClassifyBLAST -o ${nucl_fa%.fasta}.tax_guide.blastx.tab ${nucl_fa%.fasta}.tax_guide.blastx.out >/dev/null 2>&1
			taxid=$( tail -n1 ${nucl_fa%.fasta}.tax_guide.blastx.tab | cut -f2 )
			efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element Lineage >> ${nucl_fa%.fasta}.tax_guide.blastx.out
			sleep 0.4s
		fi
		if [ ! -s ${nucl_fa%.fasta}.tax_guide.blastx.out ] ; then
			echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out
		fi
		if grep -i -q "Caudovir\|Ackermannvir\|Herellevir\|Corticovir\|Levivir\|Tectivir\|crAss-like virus\|CrAssphage\|Cyanophage\|Microvir\microphage\|Siphoviridae\|Myoviridae\|phage\|Podovir\|Halovir\|sphaerolipovir\|pleolipovir\|plasmid\|Inovir\|Ampullavir\|Bicaudavir\|Fusellovir\|Guttavir\|Ligamenvir\|Plasmavir\|Salterprovir\|Cystovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			echo ${nucl_fa%.fasta}.rotate.fasta >> DTR_seqs_for_phanotate.txt
		else
			echo ${nucl_fa%.fasta}.rotate.fasta >> DTR_seqs_for_prodigal.txt
		fi
	done
	if [ -s DTR_seqs_for_phanotate.txt ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running PHANOTATE, annotate DTR contigs " $MDYT
		cat DTR_seqs_for_phanotate.txt | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/PHANOTATE/phanotate.py -f fasta -o {}.phan.fasta {}.rotate.fasta
		for PHAN in *.phan.fasta ; do 
			if [ "$ENFORCE_START_CODON" == "True" ] ; then
				sed 's/ /@/g' ${PHAN} | bioawk -c fastx '{ print }' | awk '{ if ($2 ~ /^[ATCG]TG/) { print ">"$1 ; print $2 }}' | sed 's/@/ /g' > ${PHAN%.fasta}.sort.fasta
			else
				sed 's/ /@/g' ${PHAN} | bioawk -c fastx '{ print }' | awk '{ print ">"$1 ; print $2 }' | sed 's/@/ /g' > ${PHAN%.fasta}.sort.fasta
			fi
			transeq -frame 1 -table 11 -sequence ${PHAN%.fasta}.sort.fasta -outseq ${PHAN%.phan.fasta}.trans.fasta >/dev/null 2>&1
			# put emboss transeq in directory 
			COUNTER=0 ;  
			bioawk -c fastx '{print}' ${PHAN%.phan.fasta}.trans.fasta | while read LINE ; do 
				START_BASE=$( echo $LINE | sed 's/.*START=\(.*\)\] \[.*/\1/' ) ; 
				ORF_NAME=$( echo $LINE | cut -d " " -f1 | sed 's/\(.*\)\.[0-9].*_1/\1/' ) ; 
				END_BASE=$( echo $LINE | cut -d " " -f1 | sed 's/.*\(\.[0-9].*_1\)/\1/' | sed 's/_1//g; s/\.//g' ) ; 
				ORIG_CONTIG=$( grep ">" ${PHAN%.phan.fasta}.rotate.fasta | cut -d " " -f2 ) ; 
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
				let COUNTER=COUNTER+1 ; 
				echo ">"${ORF_NAME}"_"${COUNTER} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${PHAN%.phan.fasta}.rotate.AA.fasta
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
	ROTATE_AAs=$( find * -maxdepth 0 -type f -name "${run_title}*rotate.AA.fasta" )
	if [ -n "$ROTATE_AAs" ] ; then
		for ROT in $ROTATE_AAs ; do 
			bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' $ROT > ${ROT%.fasta}.sorted.fasta
		done
	fi
fi

#-# 4 hhmscan circles/DTRs
ROTATE_SORT_AAs=$( find * -maxdepth 0 -type f -name "${run_title}*rotate.AA.sorted.fasta" )
if [ -n "$ROTATE_SORT_AAs" ] ; then
	cat $( find * -maxdepth 0 -type f -name "${run_title}*rotate.AA.sorted.fasta" ) > all_DTR_sort_genome_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_DTR_sort_genome_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_sort_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_sort_genome_proteins.AA.fasta
		SPLIT_AA_DTR_sort=$( find * -maxdepth 0 -type f -name "SPLIT_DTR_sort_GENOME_AA_*.fasta" )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running hmmscan1, annotate DTR contigs " $MDYT
		if  [[ $virus_domain_db = "standard" ]] ; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan_replicate.out --cpu 1 -E 1e-15 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_replication_clusters3 {}.fasta >/dev/null 2>&1
		elif [[ $virus_domain_db = "virion" ]]; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/virus_specific_baits_plus_missed6a {}.fasta >/dev/null 2>&1		
		elif [[ $virus_domain_db = "rna_virus" ]]; then
			echo "$SPLIT_AA_DTR_sort" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/rna_virus_rdrp_capsid_hmms1 {}.fasta >/dev/null 2>&1
		else
			echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try standard, virion, or rna_virus as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
			rm -f ./*{0..9}.fasta
			break
		fi
		HMM_REP_NUMEBR=$( find * -maxdepth 0 -type f -name "SPLIT_DTR_sort_GENOME_AA_*AA.hmmscan_replicate.out" | wc -l )
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
		if [ -s SPLIT_DTR_COMBINED.AA.hmmscan.sort.out ] ; then
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
			cp ${ROT} ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan1.fasta
		fi
	done	
	DTR_AA_FOR_HMM2=$( find * -maxdepth 0 -type f -name "${run_title}*AA.no_hmmscan1.fasta" )
	if [ -n "$DTR_AA_FOR_HMM2" ] ; then
		cat $( find * -maxdepth 0 -type f -name "${run_title}*AA.no_hmmscan1.fasta" ) > all_DTR_HMM2_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_DTR_HMM2_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_HMM2_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_HMM2_proteins.AA.fasta
			SPLIT_DTR_HMM2=$( find * -maxdepth 0 -type f -name "SPLIT_DTR_HMM2_GENOME_AA_*.fasta" )
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running hmmscan2, annotate DTR contigs " $MDYT
			echo "$SPLIT_DTR_HMM2" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan2.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/useful_hmms_baits_and_not2a {}.fasta >/dev/null 2>&1
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
	for ROT in $ROTATE_SORT_AAs ; do 
		if [ ! -s ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan2.fasta ] ; then
			cp ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan1.fasta ${ROT%.rotate.AA.sorted.fasta}.rotate.AA.no_hmmscan2.fasta
		fi
	done
	for ROT_AAs in $ROTATE_SORT_AAs ; do
		echo ">Feature "${ROT_AAs%.rotate.AA.sorted.fasta}" Table1" > ${ROT_AAs%.rotate.AA.sorted.fasta}.SCAN.tbl
		CALL_ALL_HMM=$( find * -maxdepth 0 -type f -regextype sed -regex "${ROT_AAs%.rotate.AA.sorted.fasta}\..*called_hmmscan.*txt" )
		if [ -n "$CALL_ALL_HMM" ] ; then
			cat $( find * -maxdepth 0 -type f -regextype sed -regex "${ROT_AAs%.rotate.AA.sorted.fasta}\..*called_hmmscan.*txt" ) > ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt
			if [ -s ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt ] ; then
				cat ${ROT_AAs%.rotate.AA.sorted.fasta}.all_called_hmmscans.txt | sed 's/ $//g' | while read LINE ; do 
					PROTEIN_INFO=$( grep "$LINE \[" ${ROT_AAs} ) ;  
					START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
					END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
					if grep -q "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan.sort.out ; then

						HMM_INFO=$( grep "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' ) ; 
					else
						HMM_INFO=$( grep "$LINE	" ${ROT_AAs%.AA.sorted.fasta}.AA.hmmscan2.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' )
					fi
					INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
					PROTEIN_NAME=$( echo $HMM_INFO | cut -d " " -f2- ) ; 
					echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""protein motif:$INFERENCEH" >> ${ROT_AAs%.rotate.AA.sorted.fasta}.SCAN.tbl ;
				done
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
		echo "$ROTATED_DTR_CONTIGS" | sed 's/.rotate.fasta//g' | xargs -n 1 -I {} -P $CPU blastn -query {}.rotate.fasta -db ${BLASTN_DB} -outfmt '6 std qlen slen' -max_target_seqs 100 -perc_identity 90 -num_threads 1 -out {}.blastn.out >/dev/null 2>&1
		for circle in $ROTATED_DTR_CONTIGS ; do
			if [ -s "${circle%.rotate.fasta}.blastn.out" ]; then
				python ${CENOTE_SCRIPT_DIR}/anicalc/anicalc.py -i ${circle%.rotate.fasta}.blastn.out -o ${circle%.rotate.fasta}.blastn_anicalc.out
				awk '{OFS="\t"}{FS="\t"}{ if (NR==1) {print $1, $2, $4, $5} else if ($4>=95 && $5>=85) {print $1, $2, $4, $5}}' ${circle%.rotate.fasta}.blastn_anicalc.out | head -n2 > ${circle%.rotate.fasta}.blastn_intraspecific.out
			fi
			if [ -s "${circle%.rotate.fasta}.blastn_intraspecific.out" ]; then
				INTRA_LINES=$( cat ${circle%.rotate.fasta}.blastn_intraspecific.out | wc -l | bc )	
				if [ "$INTRA_LINES" -ge 2 ] ; then
					ktClassifyBLAST -o ${circle%.rotate.fasta}.tax_guide.blastn.tab ${circle%.rotate.fasta}.blastn_intraspecific.out >/dev/null 2>&1
					taxid=$( grep -v "qname" ${circle%.rotate.fasta}.tax_guide.blastn.tab | tail -n+2 | head -n1 | cut -f2 )
					efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -tab "\n" -element Lineage ScientificName > ${circle%.rotate.fasta}.tax_guide.blastn.out
					sleep 1s
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

PROTEIN_NO_HMMSCAN2=$( find * -maxdepth 0 -type f -name "*.rotate.AA.no_hmmscan2.fasta" )

if [ -n "$PROTEIN_NO_HMMSCAN2" ]; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running RPSBLAST, annotate DTR contigs" $MDYT
	cat $( find * -maxdepth 0 -type f -name "*.rotate.AA.no_hmmscan2.fasta" ) > all_DTR_rps_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_DTR_rps_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then
	AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_DTR_RPS_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_DTR_rps_proteins.AA.fasta
		SPLIT_DTR_AA_RPS=$( find * -maxdepth 0 -type f -name "SPLIT_DTR_RPS_AA_*.fasta" )

		echo "$SPLIT_DTR_AA_RPS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t rpsblast -evalue 1e-4 -num_descriptions 5 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/cdd_rps_db/Cdd -seg yes -query {}.fasta -line_length 200 -out {}.rpsb.out >/dev/null 2>&1
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
		fi
	done
fi

# remove ORFs within ORFs that are 'hypothetical'
if [ "$ORF_WITHIN" == "True" ] ; then
	echo "$(tput setaf 5) Removing ORFS within ORFs that are 'hypothetical' $(tput sgr 0)"

	INT_TBL=$( find * -maxdepth 0 -type f -name "*.int.tbl" )
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
	INT_TBL=$( find * -maxdepth 0 -type f -name "*.int.tbl" )
	if [ -n "$INT_TBL" ] ; then
		for feat_tbl3 in $INT_TBL ; do
			sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g ; s/length=.*//g' $feat_tbl3 | sed '/--/d' > ${feat_tbl3%.int.tbl}.int2.tbl ; 
		done
	fi
fi

# Grabbing ORFs wihout RPSBLAST hits and separating them into individual files for HHsearch
echo "$(tput setaf 5) Grabbing ORFs wihout RPS-BLAST hits and separating them into individual files for HHsearch $(tput sgr 0)"

INT2_TBL=$( find * -maxdepth 0 -type f -name "*.int2.tbl" )
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
dark_orf_list=$( find * -maxdepth 0 -type f -name "*.for_hhpred.fasta" )
#-#- can this be parallelized?
if [ -n "$dark_orf_list" ] ; then
	if  [[ $HHSUITE_TOOL = "hhsearch" ]] ; then
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/hh-suite/build/src/hhsearch -i {}.for_hhpred.fasta -d $PDB_HHSUITE -d $PFAM_HHSUITE -d $CD_HHSUITE -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
			cat ${dark_orf%.for_hhpred.fasta}.out.hhr >> ${dark_orf%.rotate*.for_hhpred.fasta}.rotate.out_all.hhr ;
			rm -f ${dark_orf%.for_hhpred.fasta}.out.hhr 
		done
	elif [[ $HHSUITE_TOOL = "hhblits" ]] ; then
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/hh-suite/build/src/hhblits -i {}.for_hhpred.fasta -d $PDB_HHSUITE -d $PFAM_HHSUITE -d $CD_HHSUITE -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
			cat ${dark_orf%.for_hhpred.fasta}.out.hhr >> ${dark_orf%.rotate*.for_hhpred.fasta}.rotate.out_all.hhr ;
			rm -f ${dark_orf%.for_hhpred.fasta}.out.hhr 
		done
	else
		echo "$(tput setaf 5) Valid option for HHsuite tool (i.e. hhsearch or hhblits) was not provided. Skipping step for "$dark_orf" $(tput sgr 0)"
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	fi
	
fi

rm -f *.rotate.AA.fasta

perl ${CENOTE_SCRIPT_DIR}/hhpredreport2tbl_mt_annotation_pipe_biowulf1_gjs_edits.pl

HH_TBL=$( find * -maxdepth 0 -type f -name "*.HH.tbl" )
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
		else
			cat $feat_tbl4 | sed '/--/d; s/TPA_asm: //g; s/TPA://g' > ${feat_tbl4%.int2.tbl}.comb3.tbl
			if [ -s "${feat_tbl4%.int2.tbl}.ITR.tbl" ] ; then
				cat ${feat_tbl4%.int2.tbl}.ITR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
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
COMB3_TBL=$( find * -maxdepth 0 -type f -name "*.comb3.tbl" )
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
COMB3_TBL=$( find * -maxdepth 0 -type f -name "*.comb3.tbl" )
if [ -n "$COMB3_TBL" ] ; then
	for feat_tbl2 in *.comb3.tbl ; do 
		if grep -i -q "CRESS\|genomovir\|circovir\|bacilladnavir\|redondovir\|nanovir\|geminivir\|smacovir" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
			echo ${feat_tbl2%.comb3.tbl}" is a CRESS virus of some kind"
		else
			CONJ_COUNT=$( grep -i "virb\|type-IV\|secretion system\|conjugation\|conjugal\|transposon\| tra[a-z] \|	tra[a-z]\|trb[b-z]\|pilus" $feat_tbl2 | grep -v "TRAF\|TRAP\|protein tyrosine phosphatase\|ttRBP\|SpoU\|transport\|central" | wc -l )
			STRUCTURAL_COUNT=$( grep -i "capsid\|terminase\|portal\|baseplate\|base plate\|tail\|collar\|zot\|zonular\|minor coat\|packaging\|	virion protein" $feat_tbl2 | wc -l )
			if [[ $CONJ_COUNT -gt 0 ]] && [[ $STRUCTURAL_COUNT == 0 ]] ; then
				TAX_ORF="Conjugative Transposon"
			elif grep -i -q "zot\|zonular" $feat_tbl2 ; then
				TAX_ORF="Inoviridae"
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
			else TAX_ORF="No_suitable_orf"


			fi
			if [ "$TAX_ORF" == "Conjugative Transposon" ] ; then
				#echo "${feat_tbl2%.comb3.tbl} looks like a conjugative transposon"
				echo $TAX_ORF > ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
			elif [ "$TAX_ORF" == "Inoviridae" ] ; then
				#echo "${feat_tbl2%.comb3.tbl} looks like an Inovirus"
				echo $TAX_ORF > ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
			elif [ "$TAX_ORF" == "No_suitable_orf" ] ; then
				echo "No suitable ORF for taxonomy found for ${feat_tbl2%.comb3.tbl}, using BLASTX result."
			else

				grep -A1 "$TAX_ORF " ${feat_tbl2%.comb3.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${feat_tbl2%.comb3.tbl}.tax_orf.fasta
				if [[ $STRUCTURAL_COUNT == 1 ]] || [[ $STRUCTURAL_COUNT -gt 1 ]] ; then
					blastp -evalue 1e-2 -outfmt "6 qseqid stitle pident evalue length" -num_threads $CPU -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/blast_DBs/virus_adinto_polinton_prot_190925 -query ${feat_tbl2%.comb3.tbl}.tax_orf.fasta -out ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out >/dev/null 2>&1
				else
					blastp -evalue 1e-2 -outfmt "6 qseqid stitle pident evalue length" -num_threads $CPU -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query ${feat_tbl2%.comb3.tbl}.tax_orf.fasta -out ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out >/dev/null 2>&1
				fi
				if [ ! -s "${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out" ]; then
					echo "unclassified virus" > ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ;
				elif grep -q "virophage" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Virophage" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
				elif grep -q "adinto" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Adintovirus" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
				elif grep -i -q "polinton" ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ; then
					echo "Polinton-like virus" >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out				
				else

					ktClassifyBLAST -o ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.tab ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out >/dev/null 2>&1
					taxid=$( tail -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.tab | cut -f2 )
					efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element Lineage >> ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
					sleep 0.4s
				fi
			fi
		fi
	done
fi


cd ${base_directory}/${run_title}

LIST_OF_ITR_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "ITR_containing_contigs/*fna" )

if [ -n "$LIST_OF_ITR_DOMAIN_CONTIGS" ] ; then
	cd ITR_containing_contigs
	LIST_OF_ITR_DOMAIN_CONTIGS=$( find * -maxdepth 0 -type f -regextype sed -regex ".*.fna" )
	echo "$LIST_OF_ITR_DOMAIN_CONTIGS" | sed 's/.fna//g' | while read ITR_SEQ ; do
		if [ "$PROPHAGE" == "True" ] ;then
			cp ${ITR_SEQ}.fna no_end_contigs_with_viral_domain/${ITR_SEQ}_vs99.fna
		else
			cp ${ITR_SEQ}.fna no_end_contigs_with_viral_domain/${ITR_SEQ}.fna
		fi
	done
	cd ${base_directory}/${run_title}
else
	echo "No ITR contigs with minimum hallmark genes found."
fi

LIST_OF_VIRAL_DOMAIN_CONTIGS=$( find * -maxdepth 1 -type f -wholename "no_end_contigs_with_viral_domain/*fna" )

if [ -n "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ; then
	. ${CENOTE_SCRIPT_DIR}/annotate_linear_contigs_2.1.2.sh
else
	echo "No linear contigs with minimum hallmark genes found."
fi

cd ${base_directory}/${run_title}

CIRCULAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -wholename "DTR_contigs_with_viral_domain/*fna" )
if [ -n "$CIRCULAR_HALLMARK_CONTIGS" ] ; then
	for CIRC in $CIRCULAR_HALLMARK_CONTIGS ; do
		sed 's/ /#/g' $CIRC | bioawk -c fastx '{print ">"$name" DTR" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
	done
fi
if [ "$PROPHAGE" == "True" ] ;then
	LINEAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -regextype sed -regex ".*_vs[0-9]\{1,2\}.fna" )
	if [ -n "$LINEAR_HALLMARK_CONTIGS" ] ; then
		for LIN in $LINEAR_HALLMARK_CONTIGS ; do
			ORIGINAL_NAME=$( head -n1 ${LIN%_vs[0-9][0-9].fna}.fna | cut -d " " -f2 )
			if [ -s ${LIN%_vs[0-9][0-9].fna}.ITR.tbl ] ; then
				sed 's/ /#/g' $LIN | bioawk -v ORI="$ORIGINAL_NAME" -c fastx '{print ">"$name" "ORI" ITR" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			else
				sed 's/ /#/g' $LIN | bioawk -v ORI="$ORIGINAL_NAME" -c fastx '{print ">"$name" "ORI" no_end_feature" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
			fi
		done
	fi
else
	if [ -n "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ; then
		for LIN in $LIST_OF_VIRAL_DOMAIN_CONTIGS ; do
			sed 's/ /#/g' $LIN | bioawk -c fastx '{print ">"$name" no_end_feature" ; print $seq}' | sed 's/#/ /g' >> final_combined_virus_sequences_${run_title}.fna
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
		if grep -q "Anellovir" $tax_info ; then
			vir_name=Anelloviridae ;
		elif grep -q "Circovirus-like" $tax_info ; then
			vir_name="CRESS virus" ;
		elif grep -q "CRESS virus" $tax_info ; then
			vir_name="CRESS virus" ;
		elif grep -q "Adenovir" $tax_info ; then
			vir_name=Adenoviridae ;
		elif grep -q "Alphasatellit" $tax_info ; then
			vir_name=Alphasatellitidae ;
		elif grep -q "Ampullavir" $tax_info ; then
			vir_name=Ampullaviridae ;
		elif grep -q "Ascovir" $tax_info ; then
			vir_name=Ascoviridae ;
		elif grep -q "Asfarvir" $tax_info ; then
			vir_name=Asfarviridae ;
		elif grep -q "Bacilladnavir" $tax_info ; then
			vir_name=Bacilladnaviridae ;
		elif grep -q "Baculovir" $tax_info ; then
			vir_name=Baculoviridae ;
		elif grep -q "Bicaudavir" $tax_info ; then
			vir_name=Bicaudaviridae ;
		elif grep -q "Bidnavir" $tax_info ; then
			vir_name=Bidnaviridae ;
		elif grep -q "Ackermannvir" $tax_info ; then
			vir_name=Ackermannviridae ;
		elif grep -q "Herellevir" $tax_info ; then
			vir_name=Herelleviridae ;
		elif grep -q "Clavavir" $tax_info ; then
			vir_name=Clavaviridae ;
		elif grep -q "Adomavir" $tax_info ; then
			vir_name=Adomaviridae ;
		elif grep -q "Corticovir" $tax_info ; then
			vir_name=Corticoviridae ;
		elif grep -q "Dinodnavir" $tax_info ; then
			vir_name=Dinodnavirus ;
		elif grep -q "Autolykivir" $tax_info ; then
			vir_name=Autolykiviridae ;
		elif grep -q "Globulovir" $tax_info ; then
			vir_name=Globuloviridae ;
		elif grep -q "Pithovir" $tax_info ; then
			vir_name=Pithoviridae ;
		elif grep -q "Pandoravir" $tax_info ; then
			vir_name=Pandoravirus ;
		elif grep -q "Fusellovir" $tax_info ; then
			vir_name=Fuselloviridae ;
		elif grep -q "Guttavir" $tax_info ; then
			vir_name=Guttaviridae ;
		elif grep -q "Hepadnavir" $tax_info ; then
			vir_name=Hepadnaviridae ;
		elif grep -q "Herpesvir" $tax_info ; then
			vir_name=Herpesvirales ;
		elif grep -q "Hytrosavir" $tax_info ; then
			vir_name=Hytrosaviridae ;
		elif grep -q "Iridovir" $tax_info ; then
			vir_name=Iridoviridae ;
		elif grep -q "Lavidavir" $tax_info ; then
			vir_name=Lavidaviridae ;
		elif grep -q "Adintovir" $tax_info ; then
			vir_name=Adintovirus ;
		elif grep -q "Lipothrixvir" $tax_info ; then
			vir_name=Lipothrixviridae ;
		elif grep -q "Rudivir" $tax_info ; then
			vir_name=Rudiviridae ;
		elif grep -q "Ligamenvir" $tax_info ; then
			vir_name=Ligamenvirales ;
		elif grep -q "Marseillevir" $tax_info ; then
			vir_name=Marseilleviridae ;
		elif grep -q "Mimivir" $tax_info ; then
			vir_name=Mimiviridae ;
		elif grep -q "Nanovir" $tax_info ; then
			vir_name=Nanoviridae ;
		elif grep -q "Nimavir" $tax_info ; then
			vir_name=Nimaviridae ;
		elif grep -q "Nudivir" $tax_info ; then
			vir_name=Nudiviridae ;
		elif grep -q "Caulimovir" $tax_info ; then
			vir_name=Caulimoviridae ;
		elif grep -q "Metavir" $tax_info ; then
			vir_name=Metaviridae ;
		elif grep -q "Pseudovir" $tax_info ; then
			vir_name=Pseudoviridae ;
		elif grep -q "Retrovir" $tax_info ; then
			vir_name=Retroviridae ;
		elif grep -q "Ovalivir" $tax_info ; then
			vir_name=Ovaliviridae ;
		elif grep -q "Parvovir" $tax_info ; then
			vir_name=Parvoviridae ;
		elif grep -q "Phycodnavir" $tax_info ; then
			vir_name=Phycodnaviridae ;
		elif grep -q "Plasmavir" $tax_info ; then
			vir_name=Plasmaviridae ;
		elif grep -q "Pleolipovir" $tax_info ; then
			vir_name=Pleolipoviridae ;
		elif grep -q "Polydnavir" $tax_info ; then
			vir_name=Polydnaviridae ;
		elif grep -q "Portoglobovir" $tax_info ; then
			vir_name=Portogloboviridae ;
		elif grep -q "Poxvir" $tax_info ; then
			vir_name=Poxviridae ;
		elif grep -q "Albetovir" $tax_info ; then
			vir_name=Albetoviridae ;
		elif grep -q "Alphatetravir" $tax_info ; then
			vir_name=Alphatetraviridae ;
		elif grep -q "Alvernavir" $tax_info ; then
			vir_name=Alvernaviridae ;
		elif grep -q "Amalgavir" $tax_info ; then
			vir_name=Amalgaviridae ;
		elif grep -q "Astrovir" $tax_info ; then
			vir_name=Astroviridae ;
		elif grep -q "Aumaivir" $tax_info ; then
			vir_name=Aumaivirus ;
		elif grep -q "Avsunviroid" $tax_info ; then
			vir_name=Avsunviroidae ;
		elif grep -q "Barnavir" $tax_info ; then
			vir_name=Barnaviridae ;
		elif grep -q "Benyvir" $tax_info ; then
			vir_name=Benyviridae ;
		elif grep -q "Birnavir" $tax_info ; then
			vir_name=Birnaviridae ;
		elif grep -q "Botourmiavir" $tax_info ; then
			vir_name=Botourmiaviridae ;
		elif grep -q "Botybirnavir" $tax_info ; then
			vir_name=Botybirnavirus ;
		elif grep -q "Bromovir" $tax_info ; then
			vir_name=Bromoviridae ;
		elif grep -q "Calicivir" $tax_info ; then
			vir_name=Caliciviridae ;
		elif grep -q "Carmotetravir" $tax_info ; then
			vir_name=Carmotetraviridae ;
		elif grep -q "Chrysovir" $tax_info ; then
			vir_name=Chrysoviridae ;
		elif grep -q "Closterovir" $tax_info ; then
			vir_name=Closteroviridae ;
		elif grep -q "Cystovir" $tax_info ; then
			vir_name=Cystoviridae ;
		elif grep -q "Deltavir" $tax_info ; then
			vir_name=Deltavirus ;
		elif grep -q "Endornavir" $tax_info ; then
			vir_name=Endornaviridae ;
		elif grep -q "Flavivir" $tax_info ; then
			vir_name=Flaviviridae ;
		elif grep -q "Hepevir" $tax_info ; then
			vir_name=Hepeviridae ;
		elif grep -q "Hypovir" $tax_info ; then
			vir_name=Hypoviridae ;
		elif grep -q "Idaeovir" $tax_info ; then
			vir_name=Idaeovirus ;
		elif grep -q "Kitavir" $tax_info ; then
			vir_name=Kitaviridae ;
		elif grep -q "Levivir" $tax_info ; then
			vir_name=Leviviridae ;
		elif grep -q "Luteovir" $tax_info ; then
			vir_name=Luteoviridae ;
		elif grep -q "Matonavir" $tax_info ; then
			vir_name=Matonaviridae ;
		elif grep -q "Megabirnavir" $tax_info ; then
			vir_name=Megabirnaviridae ;
		elif grep -q "Narnavir" $tax_info ; then
			vir_name=Narnaviridae ;
		elif grep -q "Nidovir" $tax_info ; then
			vir_name=Nidovirales ;
		elif grep -q "Nodavir" $tax_info ; then
			vir_name=Nodaviridae ;
		elif grep -q "Papanivir" $tax_info ; then
			vir_name=Papanivirus ;
		elif grep -q "Partitivir" $tax_info ; then
			vir_name=Partitiviridae ;
		elif grep -q "Permutotetravir" $tax_info ; then
			vir_name=Permutotetraviridae ;
		elif grep -q "Picobirnavir" $tax_info ; then
			vir_name=Picobirnaviridae ;
		elif grep -q "Dicistrovir" $tax_info ; then
			vir_name=Dicistroviridae ;
		elif grep -q "Iflavir" $tax_info ; then
			vir_name=Iflaviridae ;
		elif grep -q "Marnavir" $tax_info ; then
			vir_name=Marnaviridae ;
		elif grep -q "Picornavir" $tax_info ; then
			vir_name=Picornaviridae ;
		elif grep -q "Polycipivir" $tax_info ; then
			vir_name=Polycipiviridae ;
		elif grep -q "Secovir" $tax_info ; then
			vir_name=Secoviridae ;
		elif grep -q "Picornavir" $tax_info ; then
			vir_name=Picornavirales ;
		elif grep -q "Pospiviroid" $tax_info ; then
			vir_name=Pospiviroidae ;
		elif grep -q "Polinton-like virus" $tax_info ; then
			vir_name="Polinton-like virus" ;
		elif grep -q "Potyvir" $tax_info ; then
			vir_name=Potyviridae ;
		elif grep -q "Quadrivir" $tax_info ; then
			vir_name=Quadriviridae ;
		elif grep -q "Reovir" $tax_info ; then
			vir_name=Reoviridae ;
		elif grep -q "Sarthrovir" $tax_info ; then
			vir_name=Sarthroviridae ;
		elif grep -q "Sinaivir" $tax_info ; then
			vir_name=Sinaivirus ;
		elif grep -q "Solemovir" $tax_info ; then
			vir_name=Solemoviridae ;
		elif grep -q "Solinvivir" $tax_info ; then
			vir_name=Solinviviridae ;
		elif grep -q "Togavir" $tax_info ; then
			vir_name=Togaviridae ;
		elif grep -q "Tombusvir" $tax_info ; then
			vir_name=Tombusviridae ;
		elif grep -q "Totivir" $tax_info ; then
			vir_name=Totiviridae ;
		elif grep -q "Tymovir" $tax_info ; then
			vir_name=Tymovirales ;
		elif grep -q "Virgavir" $tax_info ; then
			vir_name=Virgaviridae ;
		elif grep -q "Virtovir" $tax_info ; then
			vir_name=Virtovirus ;
		elif grep -q "Salterprovir" $tax_info ; then
			vir_name=Salterprovirus ;
		elif grep -q "Smacovir" $tax_info ; then
			vir_name=Smacoviridae ;
		elif grep -q "Sphaerolipovir" $tax_info ; then
			vir_name=Sphaerolipoviridae ;
		elif grep -q "Spiravir" $tax_info ; then
			vir_name=Spiraviridae ;
		elif grep -q "Crucivir" $tax_info ; then
			vir_name=Cruciviridae ;
		elif grep -q "Tectivir" $tax_info ; then
			vir_name=Tectiviridae ;
		elif grep -q "Tolecusatellit" $tax_info ; then
			vir_name=Tolecusatellitidae ;
		elif grep -q "Tristromavir" $tax_info ; then
			vir_name=Tristromaviridae ;
		elif grep -q "Turrivir" $tax_info ; then
			vir_name=Turriviridae ;
		elif grep -q "crAss-like virus\|CrAssphage" $tax_info ; then
			vir_name="CrAss-like virus" ;
		elif grep -q "Mavir\|virophage" $tax_info ; then
			vir_name=Virophage ;
		elif grep -q "Microvir" $tax_info ; then
			vir_name=Microviridae ;
		elif grep -q "microphage" $tax_info ; then
			vir_name=Microviridae ;
		elif grep -q "uncultured marine virus" $tax_info ; then
			vir_name="Virus" ;
		elif grep -q "Inovir" $tax_info ; then
			vir_name=Inoviridae ;
		elif grep -q "Siphovir" $tax_info ; then
			vir_name=Siphoviridae ;
		elif grep -q "Myovir" $tax_info ; then
			vir_name=Myoviridae ;		
		elif grep -q "unclassified dsDNA phage" $tax_info ; then
			vir_name="Phage" ;
		elif grep -q "unclassified ssDNA virus" $tax_info ; then
			vir_name="CRESS virus" ;
		elif grep -q "Lake Sarah" $tax_info ; then
			vir_name="CRESS virus" ;
		elif grep -q "Avon-Heathcote" $tax_info ; then
			vir_name="CRESS virus" ;
		elif grep -q "Circovir" $tax_info ; then
			vir_name=Circoviridae ;
		elif grep -q "Genomovir" $tax_info ; then
			vir_name=Genomoviridae ;
		elif grep -q "Geminivir" $tax_info ; then
			vir_name=Geminiviridae ;
		elif grep -q "Polyoma" $tax_info ; then
			vir_name=Polyomaviridae ;
		elif grep -q "Papillomavir" $tax_info ; then
			vir_name=Papillomaviridae ;
		elif grep -q "Halovir" $tax_info ; then
			vir_name=Halovirus ;
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
		elif grep -q "Podovir" $tax_info ; then
			vir_name=Podoviridae ;
		elif grep -q "Caudovir" $tax_info ; then
			vir_name=Caudovirales ;
		elif grep -q "dsRNA virus" $tax_info ; then
			vir_name="dsRNA virus" ;
		elif grep -q "ssRNA virus" $tax_info ; then
			vir_name="ssRNA virus" ;
		elif grep -q "unclassified RNA virus" $tax_info ; then
			vir_name="unclassified RNA viruse" ;
		elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
			vir_name="unclassified ssDNA bacterial virus" ;
		elif grep -q "phage" $tax_info ; then
			vir_name="Phage" ;
		elif grep -q "plasmid" $tax_info ; then
			vir_name="metagenomic plasmid" ;
		elif grep -q "Bacteria" $tax_info ; then
			vir_name="Phage" ;
		elif grep -q "unclassified virus" $tax_info ; then
			vir_name="virus" ;		
		elif grep -q "virus" $tax_info ; then
			vir_name="Virus" ;
		else
			if  [ -s ITR_containing_contigs/${JUST_TBL2_FILE%.comb3.tbl}.fna ] ; then
				vir_name="unclassified element" ;
			else
				vir_name="Circular genetic element" ;
			fi
		fi
		#echo $vir_name ;
		fsa_head=$( echo $vir_name " sp." )
		tax_guess=$( tail -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ) ; 
		perc_id=$( head -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out | sed 's/ /-/g' | awk '{FS="\t"; OFS="\t"} {print $2" "$3}' | sed 's/-/ /g' ) ;
		rand_id=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )

		# Editing and transferring tbl file and fasta (fsa) files to sequin directory

		if [ -s ${feat_tbl2%.comb3.tbl}.phan.fasta ]; then
			GCODE="11"
		else
			GCODE="1"
		fi
		if echo "$feat_tbl2" | grep -q "DTR_contigs_with_viral_domain" ; then
			TOPOLOGY="circular"
			NUCL_FILE="${feat_tbl2%.comb3.tbl}.rotate.fasta"
			input_contig_name=$( head -n1 ${feat_tbl2%.comb3.tbl}.rotate.fasta | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' ) 
		else
			TOPOLOGY="linear"
			NUCL_FILE="${feat_tbl2%.comb3.tbl}.fna"
			if [ "$PROPHAGE" == "True" ] ; then
				input_contig_name=$( head -n1 ${feat_tbl2%_vs[0-9][0-9].comb3.tbl}.fna | cut -d " " -f 2 | sed 's/|.*//g; s/>//g' )
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
		bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -v topoq="$TOPOLOGY" -v gcodeq="$GCODE" -v o_name="$input_contig_name" -v crispr1="$CRISPR" -v blastn="$BLASTN_INFO" -c fastx '{ print ">" newname " [note=input name:"o_name" -- closest relative: " tax_var " " perc_var " ; " crispr1" "blastn"] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology="topoq"] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode="gcodeq"]" ; print $seq }' $NUCL_FILE > sequin_and_genome_maps/${JUST_TBL2_FILE%.comb3.tbl}.fsa ; 

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
	COMB4_TBL=$( find * -maxdepth 0 -type f -name "*.tbl" )
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
	cd ..
fi
### conjugative machinery table
if [ -d sequin_and_genome_maps ] ; then
	cd sequin_and_genome_maps
	COMB4_TBL=$( find * -maxdepth 0 -type f -name "*.tbl" )
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
	cd ..
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
CIRCULAR_HALLMARK_CONTIGS=$( find * -maxdepth 1 -type f -wholename "DTR_contigs_with_viral_domain/*fna" )

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
			ORIGINAL_NAME=$( head -n1 ${LINEAR%_vs[0-9][0-9].fna}.fna | cut -d " " -f2 )
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
		END_FEATURE="None"
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
	cd ..
fi
rm -rf bt2_indices/
rm -f other_contigs/*.AA.fasta other_contigs/*.AA.sorted.fasta other_contigs/*.out other_contigs/*.dat other_contigs/*called_hmmscan.txt other_contigs/SPLIT_LARGE_GENOME_AA_*fasta ITR_containing_contigs/SPLIT_ITR_AA*fasta SPLIT_CIRCULAR_AA* *called_hmmscan.txt circular_contigs_spades_names.txt
rm -f no_end_contigs_with_viral_domain/*.called_hmmscan2.txt no_end_contigs_with_viral_domain/*.hmmscan2.out no_end_contigs_with_viral_domain/*all_hhpred_queries.AA.fasta no_end_contigs_with_viral_domain/*.all_start_stop.txt no_end_contigs_with_viral_domain/*.trnascan-se2.txt no_end_contigs_with_viral_domain/*.for_hhpred.txt no_end_contigs_with_viral_domain/*.for_blastp.txt no_end_contigs_with_viral_domain/*.HH.tbl no_end_contigs_with_viral_domain/*.hypo_start_stop.txt  no_end_contigs_with_viral_domain/*.remove_hypo.txt no_end_contigs_with_viral_domain/*.rps_nohits.fasta no_end_contigs_with_viral_domain/*.tax_guide.blastx.tab no_end_contigs_with_viral_domain/*.tax_orf.fasta no_end_contigs_with_viral_domain/*.trans.fasta no_end_contigs_with_viral_domain/*.called_hmmscan*.txt no_end_contigs_with_viral_domain/*.no_hmmscan*.fasta  no_end_contigs_with_viral_domain/SPLIT_LIN_HMM2_GENOME_AA*fasta no_end_contigs_with_viral_domain/SPLIT_LIN_sort_GENOME_AA* no_end_contigs_with_viral_domain/SPLIT_LIN_RPS_AA* no_end_contigs_with_viral_domain/*used_positions.txt no_end_contigs_with_viral_domain/*seq_chunk_coordinates.csv no_end_contigs_with_viral_domain/*blast_hypo.fasta no_end_contigs_with_viral_domain/*CDS.bed no_end_contigs_with_viral_domain/*tRNA.bed no_end_contigs_with_viral_domain/*ORFs_over_tRNAs.tsv no_end_contigs_with_viral_domain/*prodigal.fasta no_end_contigs_with_viral_domain/*all_called_hmmscans.txt no_end_contigs_with_viral_domain/*phan.fasta no_end_contigs_with_viral_domain/*phan.sort.fasta

echo "$(tput setaf 3)output directory: "$run_title" $(tput sgr 0)"
echo "$(tput setaf 3) >>>>>>CENOTE-TAKER 2 HAS FINISHED TAKING CENOTES<<<<<< $(tput sgr 0)"



