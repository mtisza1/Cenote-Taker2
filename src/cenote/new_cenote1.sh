#!/bin/bash

## Setting input parameters  

original_contigs=$1
run_title=$2

base_directory=$PWD

circ_length_cutoff=1000
linear_length_cutoff=1000
virus_domain_db="virion"
LIN_MINIMUM_DOMAINS=2
CIRC_MINIMUM_DOMAINS=1
MEM=25
CPU=12
WRAP="True"

CENOTE_DBs="/Users/u241374/mike_tisza/cmmr_repos/Cenote-Taker2"
CENOTE_SCRIPTS="/Users/u241374/mike_tisza/cmmr_repos/Cenote-Taker2/src/cenote"

echo "@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "Your specified arguments:"
echo "original contigs:                  $original_contigs"
echo "title of this run:                 $run_title"
echo "minimum circular contig length:    $circ_length_cutoff"
echo "minimum linear contig length:      $linear_length_cutoff"
echo "virus domain database:             $virus_domain_db"
echo "min. viral hallmarks for linear:   $LIN_MINIMUM_DOMAINS"
echo "min. viral hallmarks for circular: $CIRC_MINIMUM_DOMAINS"
echo "GB of memory:                      $MEM"
echo "number of CPUs available for run:  $CPU"
echo "Cenote DBs directory:              $CENOTE_DBs"
echo "Cenote scripts directory:          $CENOTE_SCRIPTS"


echo " "

MDYT=$( date +"%m-%d-%y---%T" )
echo "time update: configuring run directory " $MDYT


#checking validity of run_title
if [[ "$run_title" =~ ^[a-zA-Z0-9_]+$ ]] && [ ${#run_title} -le 18 ] ; then 
	echo $run_title ; 
else
	echo "$run_title is not a valid name for the run title ( -r argument)"
	echo " the run title needs to be only letters, numbers and underscores (_) and 18 characters or less. Exiting."
	exit
fi

if [ -s ${original_contigs} ] ; then 
	original_con_base=$( basename $original_contigs )  
else  
	echo "${original_contigs} not found"
	exit
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

if [ ! -d "${run_title}/ct2_tmp" ]; then
	mkdir ${run_title}/ct2_tmp
fi

TEMP_DIR="${run_title}/ct2_tmp"

# Removing contigs under $circ_length_cutoff nts and detecting circular contigs
if [ $circ_length_cutoff -gt $linear_length_cutoff ] ; then
	LENGTH_MINIMUM=$linear_length_cutoff
else
	LENGTH_MINIMUM=$circ_length_cutoff

fi

if [ -s ${original_contigs} ] ; then
	original_con_base=$( basename $original_contigs )

	## filter on minimum length
	bioawk -v run_var="$run_title" -v contig_cutoff="$LENGTH_MINIMUM" -c fastx\
	  '{ if(length($seq) > contig_cutoff) { print ">"run_var NR" "$name; print $seq }}' \
	  $original_contigs > ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta


else
	echo "${original_contigs} not found"
	exit
fi


## split contigs into equal parts for prodigal ORF calling
if [ -s ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta ] ; then
	if [ ! -d "${TEMP_DIR}/split_orig_contigs" ]; then
		mkdir ${TEMP_DIR}/split_orig_contigs
	fi

	seqkit split --quiet -j $CPU -p $CPU -O ${TEMP_DIR}/split_orig_contigs ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta

	SPLIT_ORIG_CONTIGS=$( find ${TEMP_DIR}/split_orig_contigs -type f -name "*.fasta" )

	if [ ! -z "$SPLIT_ORIG_CONTIGS" ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running prodigal on all contigs " $MDYT

		echo "$SPLIT_ORIG_CONTIGS" | sed 's/.fasta//g' |\
		  xargs -n 1 -I {} -P $CPU -t prodigal -a {}.prod.faa -i {}.fasta -p meta -q >/dev/null 2>&1

	else
		echo "can't find split original contigs"

	fi
else
	echo "couldn't find ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta"
	echo "exiting"
	exit
fi


## run pyhmmer on prodigal ORF files
SPLIT_ORIG_AAs=$( find ${TEMP_DIR}/split_orig_contigs -type f -name "*.prod.faa" )

if [ ! -z "$SPLIT_ORIG_CONTIGS" ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running pyhmmer on all ORFs " $MDYT

	python ${CENOTE_SCRIPTS}/python_modules/pyhmmer_runner.py ${TEMP_DIR}/split_orig_contigs ${TEMP_DIR}/split_orig_pyhmmer\
	  ${CENOTE_DBs}/hmmscan_DBs/virus_specific_baits_plus_missed6a.h3m
else
	echo "couldn't find prodigal AA seqs in ${TEMP_DIR}/split_orig_contigs"
fi

## keep contigs with minimum hallmark genes

# set minimum hallmark genes
if [ $CIRC_MINIMUM_DOMAINS -gt $LIN_MINIMUM_DOMAINS ] ; then
	HALLMARK_MINIMUM=$LIN_MINIMUM_DOMAINS
else
	HALLMARK_MINIMUM=$CIRC_MINIMUM_DOMAINS
fi

## grabbing contigs with minimum marker gene number
if [ -s ${TEMP_DIR}/split_orig_pyhmmer/contig_hallmark_count.tsv ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "grabbing contigs with minimum marker gene number " $MDYT

	KEEP_COUNT=$( tail -n+2 ${TEMP_DIR}/split_orig_pyhmmer/contig_hallmark_count.tsv | awk -v min="${HALLMARK_MINIMUM}" '{OFS=FS="\t"}{if ($2 >= min) {print}}' | wc -l )

	if [ $KEEP_COUNT -ge 1 ] ; then
		echo "found ${KEEP_COUNT} contig(s) with at least ${HALLMARK_MINIMUM} hallmark gene(s)"
		echo "analyzing these contigs further..."

		tail -n+2 ${TEMP_DIR}/split_orig_pyhmmer/contig_hallmark_count.tsv |\
		  awk -v min="${HALLMARK_MINIMUM}" '{OFS=FS="\t"}{if ($2 >= min) {print $1}}' > ${TEMP_DIR}/split_orig_pyhmmer/contigs_to_keep.txt

		seqkit grep -f ${TEMP_DIR}/split_orig_pyhmmer/contigs_to_keep.txt\
		  ${run_title}/${run_title}.contigs_over_${LENGTH_MINIMUM}nt.fasta > ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta

	else
		echo "no contigs with ${KEEP_COUNT} hallmark gene(s) were detected. Exiting."
		exit
	fi
else
	echo "couldn't find ${TEMP_DIR}/split_orig_pyhmmer/contig_hallmark_count.tsv. Exiting."
	exit
fi

## detecting DTRs and ITRs. Trimming DTRs
if [ -s ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta ] ; then

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "detecting DTRs and ITRs " $MDYT

	python ${CENOTE_SCRIPTS}/python_modules/terminal_repeats.py ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta ${TEMP_DIR} $WRAP

else
	echo "couldn't find hallmark contigs at ${TEMP_DIR}/unprocessed_hallmark_contigs.fasta"
	exit
fi

## Rotating DTR contigs
if [ -s ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta ] && [ -s ${TEMP_DIR}/hallmark_contigs_terminal_repeat_summary.tsv ] ; then
	if [ "$WRAP" == "True" ] ; then

		MDYT=$( date +"%m-%d-%y---%T" )
		echo "rotating DTR contigs " $MDYT

		if [ ! -d "${TEMP_DIR}/rotation" ]; then
			mkdir ${TEMP_DIR}/rotation
		fi

		awk '{OFS=FS="\t"}{ if (NR != 1 && $4 != "NA") {print $1, $3} }' ${TEMP_DIR}/hallmark_contigs_terminal_repeat_summary.tsv |\
		  while read DTR_CONTIG TRIM_LENGTH ; do
			seqkit grep -p "${DTR_CONTIG}" ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta |\
			  prodigal -a ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa -i /dev/stdin -c -p meta -q >/dev/null 2>&1
			FWD_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
				awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
			REV_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
				awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == -1) {print $4}}' | wc -l )
			if [ $FWD_GENES -ge $REV_GENES ] && [ $FWD_GENES -ge 1 ]; then
				START_BASE=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.unrotated.faa | sed 's/ # /	/g' |\
					awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' |\
					sort -rg -k2,2 | head -n1 | cut -f1 )
				seqkit grep -p "${DTR_CONTIG}" ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta |\
				  seqkit restart -i ${START_BASE} > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
			elif [ $REV_GENES -ge 1 ]; then
				seqkit grep -p "${DTR_CONTIG}" ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta |\
				  seqkit seq --quiet -t DNA -r -p > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna
				prodigal -a ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa -i ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna -p meta -q >/dev/null 2>&1
				RC_FWD_GENES=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa | sed 's/ # /	/g' |\
					awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $4}}' | wc -l )
				if [ $RC_FWD_GENES -ge 1 ] ; then 
					START_BASE=$( grep "^>" ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.faa | sed 's/ # /	/g' |\
						awk '{FS=OFS="\t"}{if ($0 ~ "partial=00;start_type" && $4 == 1) {print $2, ($3-$2)}}' |\
						sort -rg -k2,2 | head -n1 | cut -f1 )
					cat ${TEMP_DIR}/rotation/${DTR_CONTIG}.rc.fna |\
					  seqkit restart -i ${START_BASE} > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
				else
					echo "Can't find suitable ORF to set rotation of ${DTR_CONTIG} and will remain unrotated"
					seqkit grep -p "${DTR_CONTIG}" ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
				fi
			else
				echo "Can't find suitable ORF to set rotation of $nucl_fa and will remain unrotated"
				seqkit grep -p "${DTR_CONTIG}" ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta > ${TEMP_DIR}/rotation/${DTR_CONTIG}.rotate.fasta
			fi
		done

		awk '{OFS=FS="\t"}{ if (NR != 1 && $4 == "NA") {print $1} }' ${TEMP_DIR}/hallmark_contigs_terminal_repeat_summary.tsv > ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt

		if [ -s ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt ] ; then
			seqkit grep -f ${TEMP_DIR}/hallmark_contigs_no_DTRS.txt\
			  ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta > ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta
		fi

		ALL_DTRS=$( find ${TEMP_DIR}/rotation -type f -name "*.rotate.fasta" )

		if [ -n "$ALL_DTRS" ] ; then
			echo "$ALL_DTRS" | while read SEQ ; do
				cat $SEQ
			done > ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta
		fi

		if [ -s ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta ] ; then
			cat ${TEMP_DIR}/hallmark_contigs_with_DTRS.rotated.fasta >> ${TEMP_DIR}/oriented_hallmark_contigs.fasta
		fi

		if [ -s ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta ] ; then
			cat ${TEMP_DIR}/hallmark_contigs_no_DTRS.fasta >> ${TEMP_DIR}/oriented_hallmark_contigs.fasta
		fi


	else
		echo "not wrapping"
		
		cp ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta ${TEMP_DIR}/oriented_hallmark_contigs.fasta
	fi
else
	echo "couldn't find hallmark contigs with processed terminal repeats at ${TEMP_DIR}/trimmed_TRs_hallmark_contigs.fasta"
fi






