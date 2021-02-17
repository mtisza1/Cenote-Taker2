#!/bin/bash

# Split from anotation_of_viral_domain_contigs.sh on 190919

cd no_end_contigs_with_viral_domain/

echo "pruning script opened"
vd_fastas=$( find * -maxdepth 0 -type f -name "*.fna" )

if [ -n "$vd_fastas" ] ; then
	echo "fna files found"
	for NO_END in $vd_fastas ; do
		LENGTH_SEQ=$( bioawk -c fastx '{print length($seq)}' $NO_END )
		if [[ "$LENGTH_SEQ" -lt 10000 ]] ; then 
			bioawk -c fastx '{print ">"$name"_vs01 "$4 ; print $seq}' $NO_END > ${NO_END%.fna}_vs01.fna
			mv ${NO_END%.fna}.AA.sorted.fasta ${NO_END%.fna}_vs01.AA.sorted.fasta
			echo "$NO_END is too short to prune chromosomal regions"
		elif [[ "$LENGTH_SEQ" -gt 10000 ]] || [[ "$LENGTH_SEQ" == 10000 ]] ; then
			mv ${NO_END%.fna}.AA.sorted.fasta ${NO_END%.fna}.AA.sorted1.fasta
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "Making table files of viral bait outputs" $MDYT

			#grep -v "^#\|plasmid_cluster" ${NO_END%.fna}.AA.hmmscan.out | sed 's/ \+/	/g' | sort -u -k3,3 > ${NO_END%.fna}.AA.hmmscan.sort.out
			if [ ! -z "${NO_END%.fna}.AA.hmmscan.sort.out" ] ; then
				cut -f3 ${NO_END%.fna}.AA.hmmscan.sort.out | awk '{ print $0" " }' > ${NO_END%.fna}.AA.called_hmmscan1.txt ; 
				cat ${NO_END%.fna}.AA.called_hmmscan1.txt | while read LINE ; do 
					PROTEIN_INFO=$( grep "$LINE \[" ${NO_END%.fna}.AA.sorted1.fasta ) ;  
					START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
					END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
					HMM_INFO=$( grep "$LINE	" ${NO_END%.fna}.AA.hmmscan.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' ) ; 
					INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
					echo -e "$LINE\t$START_BASEH\t$END_BASEH\t$INFERENCEH\t$HMM_INFO"
				done > ${NO_END%.fna}.VIRUS_BAIT_TABLE.txt
			fi

		fi
	done


	### redo hmmscan parallelization
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: HMMSCAN of common viral domains beginning" $MDYT
	cat $( find * -maxdepth 0 -type f -name "*.AA.sorted1.fasta" ) > all_prunable_seq_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_prunable_seq_proteins.AA.fasta | wc -l | bc )
	AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
	if [ $AA_SEQS_PER_FILE = 0 ] ; then
		AA_SEQS_PER_FILE=1
	fi
	awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_PRUNE_SEQ_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_prunable_seq_proteins.AA.fasta
	SPLIT_AA_PRUNE=$( find * -maxdepth 0 -type f -name "SPLIT_PRUNE_SEQ_AA_*.fasta" )
	echo "$SPLIT_AA_PRUNE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan2.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/useful_hmms_baits_and_not2a {}.fasta >/dev/null 2>&1
	cat SPLIT_PRUNE_SEQ_AA_*AA.hmmscan2.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_PRUNE_SEQ_COMBINED.AA.hmmscan2.sort.out
	if [ -s SPLIT_PRUNE_SEQ_COMBINED.AA.hmmscan2.sort.out ] ; then
		cut -f3 SPLIT_PRUNE_SEQ_COMBINED.AA.hmmscan2.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
			### not sure on this
			grep "${HIT}_" SPLIT_PRUNE_SEQ_COMBINED.AA.hmmscan2.sort.out > ${HIT}.AA.hmmscan2.sort.out
			# ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.no_hmmscan1.fasta
			grep "${HIT}_" SPLIT_PRUNE_SEQ_COMBINED.AA.hmmscan2.sort.out | sort -u -k3,3 | cut -f3 > ${HIT}.AA.called_hmmscan2.txt
			grep -v -f ${HIT}.AA.called_hmmscan2.txt ${HIT}.AA.sorted1.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.AA.prune_no_hmmscan2.fasta

		done
	fi
	CALLED_VIRAL=$( find * -maxdepth 0 -type f -name "*.AA.called_hmmscan2.txt" )
	if [ -n "$CALLED_VIRAL" ] ; then
		for NO_END in $CALLED_VIRAL ; do
			cat $NO_END | while read LINE ; do 
				PROTEIN_INFO=$( grep "$LINE \[" ${NO_END%.AA.called_hmmscan2.txt}.AA.sorted1.fasta ) ;  
				START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
				END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
				HMM_INFO=$( grep "$LINE	" ${NO_END%.AA.called_hmmscan2.txt}.AA.hmmscan2.sort.out | head -n1 | cut -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g' ) ; 
				INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
				PROTEIN_NAME=$( echo $HMM_INFO | cut -d " " -f2- ) ; 
				if [[ $INFERENCEH == *"UniProtKB"* ]]; then
					echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""similar to AA sequence:$INFERENCEH" >> ${NO_END%.AA.called_hmmscan2.txt}.SCAN.tbl ; 
				else
					echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""protein motif:$INFERENCEH" >> ${NO_END%.AA.called_hmmscan2.txt}.SCAN.tbl ; 
				fi
			done
		done
	fi
	SCAN1_TBL=$( find * -maxdepth 0 -type f -name "*.SCAN.tbl" )
	if [ -n "$SCAN1_TBL" ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: making tables for hmmscan and rpsblast outputs " $MDYT
		for feat_tbl2 in $SCAN1_TBL ; do
			SCAN_TBL_LENGTH=$( cat $feat_tbl2 | wc -l )
			if [[ "$SCAN_TBL_LENGTH" -gt 1 ]] ; then
				echo "making scan table for $feat_tbl2"
				grep "^[0-9]" -A3 $feat_tbl2 | sed '/--/d' | sed 's/ /_/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n		//g' | while read LINE ; do
					if echo $LINE | grep -q "CDS" ; then
						GENOME=${feat_tbl2%.SCAN.tbl}
						FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
						FEAT_START=$( echo $LINE | cut -d " " -f1 )
						FEAT_END=$( echo $LINE | cut -d " " -f2 )
						FEAT_NAME=$( echo $LINE | cut -d " " -f7 )
						FEAT_ATT=$( echo $LINE | cut -d " " -f9 )
						FEAT_ID=$( echo $LINE | cut -d " " -f5 | sed 's/lcl|//g' )
						echo -e "$FEAT_ID\t$FEAT_START\t$FEAT_END\tVIRUS_COMMON\t$FEAT_NAME"
					fi
				done > ${feat_tbl2%.SCAN.tbl}.HMMSCAN_TABLE.txt
			else
				echo "not making scan table for $feat_tbl2"
			fi
		done
	fi
	### redo RPS part
	NO_HMMSCAN_AA=$( find * -maxdepth 0 -type f -name "*.AA.prune_no_hmmscan2.fasta" )
	if [ -n "$NO_HMMSCAN_AA" ] ; then
		cat $( find * -maxdepth 0 -type f -name "*.AA.prune_no_hmmscan2.fasta" ) > all_prunable_rps_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_prunable_rps_proteins.AA.fasta | wc -l | bc )
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_PRUNE_RPS_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_prunable_rps_proteins.AA.fasta
		SPLIT_AA_RPS=$( find * -maxdepth 0 -type f -name "SPLIT_PRUNE_RPS_AA_*.fasta" )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running RPSBLAST on each sequence " $MDYT
		echo "$SPLIT_AA_RPS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t rpsblast -evalue 1e-4 -num_descriptions 5 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/cdd_rps_db/Cdd -seg yes -query {}.fasta -line_length 200 -out {}.rpsb.out >/dev/null 2>&1
		cat *rpsb.out > COMBINED_RESULTS_PRUNE.AA.rpsblast.out
		perl ${CENOTE_SCRIPT_DIR}/rpsblastreport_to_table2.pl
		cut -f1 COMBINED_RESULTS_PRUNE.RPS_TABLE.txt | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read CONTIG ; do
			grep "^${CONTIG}_" COMBINED_RESULTS_PRUNE.RPS_TABLE.txt > ${CONTIG}.RPS_TABLE.txt
		done
	fi

	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: parsing tables into virus_signal.seq files for hmmscan and rpsblast outputs " $MDYT
	HMM_TBL=$( find * -maxdepth 0 -type f -name "*.HMMSCAN_TABLE.txt" )
	if [ -n "$HMM_TBL" ] ; then
		for TABLE1 in $HMM_TBL ; do
			echo $TABLE1
			CONTIG_LENGTH=$( bioawk -c fastx '{ print length($seq) }' ${TABLE1%.HMMSCAN_TABLE.txt}.fna )
			awk -v lengthq="$CONTIG_LENGTH" 'BEGIN {for (i=1;i<=lengthq;i++) print i, "Z"}' | sed 's/ /	/g' > ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab


			cat ${TABLE1%.HMMSCAN_TABLE.txt}.RPS_TABLE.txt | while read LINE ; do
				RPS_START=$( echo "$LINE" | cut -f2 ) ;
				RPS_END=$( echo "$LINE" | cut -f3 ) ; 
				INFER=$( echo "$LINE" | cut -f4 ) ;
				if [ $INFER == "HYPO" ] ; then
					if [[ "$RPS_END" -gt "$RPS_START" ]] ; then
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=startq;i<=endq;i++) print i, "X"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					elif [[ "$RPS_START" -gt "$RPS_END" ]] ; then 
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=endq;i<=startq;i++) print i, "X"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					fi
				elif [[ $INFER == *"PHA0"* ]] || grep -q "$INFER" ${CENOTE_SCRIPT_DIR}/viral_cdds_and_pfams_191028.txt ; then
					if [[ "$RPS_END" -gt "$RPS_START" ]] ; then 
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=startq;i<=endq;i++) print i, "V"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					elif [[ "$RPS_START" -gt "$RPS_END" ]] ; then 
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=endq;i<=startq;i++) print i, "V"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					fi	
				else
					if [[ "$RPS_END" -gt "$RPS_START" ]] ; then 
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=startq;i<=endq;i++) print i, "Y"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					elif [[ "$RPS_START" -gt "$RPS_END" ]] ; then 
						awk -v startq="$RPS_START" -v endq="$RPS_END" 'BEGIN {for (i=endq;i<=startq;i++) print i, "Y"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab
					fi
				fi
			done	

			cat $TABLE1 | while read LINE ; do 
				VIR_START=$( echo "$LINE" | cut -f2 ) ; 
				VIR_END=$( echo "$LINE" | cut -f3 ) ; 
				if [[ "$VIR_END" -gt "$VIR_START" ]] ; then 
					awk -v startq="$VIR_START" -v endq="$VIR_END" 'BEGIN {for (i=startq;i<=endq;i++) print i, "V"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab

				elif [[ "$VIR_START" -gt "$VIR_END" ]] ; then 
					awk -v startq="$VIR_START" -v endq="$VIR_END" 'BEGIN {for (i=endq;i<=startq;i++) print i, "V"}' | sed 's/ /	/g' >> ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab 
				fi ; 
			done

			sort -g -k1,1 ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.tab | awk '!_[$1]++' | cut -f2 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n//g' > ${TABLE1%.HMMSCAN_TABLE.txt}.virus_signal.seq

		done
	fi
else
	echo "fna files NOT found"
fi
# Run Anna's python for extracting virus segments from potential prophage
if [ $PROPHAGE == "False" ] ; then
	echo "prophages not being pruned"

else
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: Identifying virus chunks, chromosomal junctions, and pruning contigs as necessary " $MDYT
	SIGNAL_SEQ=$( find * -maxdepth 0 -type f -name "*virus_signal.seq" )
	if [ -n "$SIGNAL_SEQ" ] ; then
		for VSEQ in *virus_signal.seq ; do
			python ${CENOTE_SCRIPT_DIR}/cenote_virus_segments_V6.py $VSEQ
		done
	else
		echo "No prophage signal files found"
	fi

	# reformat .fastas to get viral segments of prophage
	CHUNK_COORDINATES=$( find * -maxdepth 0 -type f -name "*.virus_signal.seq_chunk_coordinates.csv" )
	if [ -n "$CHUNK_COORDINATES" ] ; then
		VIR_COUNTER=0 
		for CHUNKS in *.virus_signal.seq_chunk_coordinates.csv ; do
			cat $CHUNKS | while read LINE ; do
				CHUNK_START=$( echo "$LINE" | cut -d "," -f2 )
				CHUNK_END=$( echo "$LINE" | cut -d "," -f3 )
				CHUNK_VIR=$( cat ${CHUNKS%.virus_signal.seq_chunk_coordinates.csv}.VIRUS_BAIT_TABLE.txt | while read BAIT_LINE ; do echo "$BAIT_LINE" | awk -v CHUNK_STARTQ="$CHUNK_START" -v CHUNK_ENDQ="$CHUNK_END" '{ if ($3-(($3-$2)/2) >= CHUNK_STARTQ && $3-(($3-$2)/2) <= CHUNK_ENDQ) { print } }' ; done | wc -l ) ;
				
				if [[ $CHUNK_VIR -gt $LIN_MINIMUM_DOMAINS ]] || [ $CHUNK_VIR == $LIN_MINIMUM_DOMAINS ] ; then 
					let VIR_COUNTER=VIR_COUNTER+1 ;
					VIR_VALUE=$( printf "%02d" $VIR_COUNTER)
					CHUNK_LENGTH=$(( ${CHUNK_END}-${CHUNK_START} ))
					bioawk -v chunk_startq="$CHUNK_START" -v chunk_endq="$CHUNK_END" -v chunk_lengthq="$CHUNK_LENGTH" -v parent="${CHUNKS%.virus_signal.seq_chunk_coordinates.csv}" -v good_chunk="$VIR_VALUE" -c fastx '{ print ">"parent"_vs"good_chunk" "chunk_startq"-"chunk_endq ; print substr($seq, chunk_startq, chunk_lengthq)}' ${CHUNKS%.virus_signal.seq_chunk_coordinates.csv}.fna > ${CHUNKS%.virus_signal.seq_chunk_coordinates.csv}_vs${VIR_VALUE}.fna
				fi
			done
		done
	fi
fi

POST_PRUNE_CONTIGS=$( find * -maxdepth 1 -type f -regextype sed -regex ".*_vs[0-9]\{1,2\}.fna" )

if [ -n "$POST_PRUNE_CONTIGS" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: Making prophage table " $MDYT
	echo "SEQUENCE_FILE	CENOTE_PARENT_NAME	INPUT_PARENT_NAME	FINAL_LENGTH	PARENT_LENGTH	PRUNING_TRIED	CHROM_REMOVED	LEFT_JUNCTION	RIGHT_JUNCTION	HALLMARK_COUNT	HALLMARK_GENES	HALLMARK_AA_NAMES" > ${run_title}_PRUNING_INFO_TABLE.tsv
	for CONTIG in $POST_PRUNE_CONTIGS ; do
		CENOTE_PARENT=${CONTIG%_vs[0-9][0-9].fna}.fna
		PARENT_LENGTH=$( bioawk -c fastx '{print length($seq)}' $CENOTE_PARENT )
		PRUNED_LENGTH=$( bioawk -c fastx '{print length($seq)}' $CONTIG )		
		if [ -s ${CONTIG%_vs[0-9][0-9].fna}.VIRUS_BAIT_TABLE.txt ] ; then
			LEFT_COORD=$( head -n1 $CONTIG | cut -d " " -f2 | cut -d "-" -f1 )
			RIGHT_COORD=$( head -n1 $CONTIG | cut -d " " -f2 | cut -d "-" -f2 )
			PRUNING_TRIED="True"
			PQ_LENGTH=$(( ${PARENT_LENGTH}+1 ))
			if [ $LEFT_COORD == 0 ] && [ $RIGHT_COORD == $PQ_LENGTH ] ; then
				CHROM_REMOVED="False"
			else
				CHROM_REMOVED="True"
			fi
			HALLMARK_COUNT=$( awk -v LEFTQ="$LEFT_COORD" -v RIGHTQ="$RIGHT_COORD" '{FS="\t"}{OFS="\t"}{ if ($2>LEFTQ && $2<RIGHTQ && $3>LEFTQ && $3<RIGHTQ) {print $5}}' ${CONTIG%_vs[0-9][0-9].fna}.VIRUS_BAIT_TABLE.txt | wc -l | bc )
			HALLMARK_GENES=$( awk -v LEFTQ="$LEFT_COORD" -v RIGHTQ="$RIGHT_COORD" '{FS="\t"}{OFS="\t"}{ if ($2>LEFTQ && $2<RIGHTQ && $3>LEFTQ && $3<RIGHTQ) {print $5}}' ${CONTIG%_vs[0-9][0-9].fna}.VIRUS_BAIT_TABLE.txt | sed 's/ /_/g ; s/,//g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' )
			HALLMARK_AA_NAMES=$( awk -v LEFTQ="$LEFT_COORD" -v RIGHTQ="$RIGHT_COORD" '{ if ($2>LEFTQ && $2<RIGHTQ && $3>LEFTQ && $3<RIGHTQ) {print $1}}' ${CONTIG%_vs[0-9][0-9].fna}.VIRUS_BAIT_TABLE.txt | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' )
		else
			LEFT_COORD="None"
			RIGHT_COORD="None"
			PRUNING_TRIED="False"
			CHROM_REMOVED="False"
			HALLMARK_COUNT=$( cat ${CONTIG%_vs[0-9][0-9].fna}.AA.hmmscan.sort.out | wc -l | bc )
			HALLMARK_GENES=$( cut -f1 ${CONTIG%_vs[0-9][0-9].fna}.AA.hmmscan.sort.out | sed 's/ /_/g ; s/,//g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' )
			HALLMARK_AA_NAMES=$( cut -f3 ${CONTIG%_vs[0-9][0-9].fna}.AA.hmmscan.sort.out | sed 's/ /_/g ; s/,//g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/\|/g' )
		fi
		ORIGINAL_NAME=$( head -n1 ${CONTIG%_vs[0-9][0-9].fna}.fna | cut -d " " -f2 )
		echo "${CONTIG}	${CENOTE_PARENT}	${ORIGINAL_NAME}	${PRUNED_LENGTH}	${PARENT_LENGTH}	${PRUNING_TRIED}	${CHROM_REMOVED}	${LEFT_COORD}	${RIGHT_COORD}	${HALLMARK_COUNT}	${HALLMARK_GENES}	${HALLMARK_AA_NAMES}" >> ${run_title}_PRUNING_INFO_TABLE.tsv
		### finish this


	done
fi

if [ -s ${run_title}_PRUNING_INFO_TABLE.tsv ] ; then
	mv ${run_title}_PRUNING_INFO_TABLE.tsv ../
fi

rm -f *.virus_signal.tab *.used_positions.txt *.phan.fasta *.phan.sort.fasta *rpsb.out SPLIT_PRUNE_RPS_AA*fasta SPLIT_PRUNE_SEQ_AA*fasta

cd ..

echo "$(tput setaf 3) FINISHED PRUNING CONTIGS WITH AT LEAST $LIN_MINIMUM_DOMAINS VIRAL DOMAIN(S) $(tput sgr 0)"
