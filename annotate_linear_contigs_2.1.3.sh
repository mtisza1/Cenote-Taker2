#!/bin/bash

### annotate linear sequences
cd no_end_contigs_with_viral_domain/

rm -f *no_hmmscan1.fasta
#CIRCULAR_HALLMARK_CONTIGS=$( find . -maxdepth 1 -type f -name "*fna" )
echo "Annotating linear contigs"

#-# blastx for translation decision
if [ "$PROPHAGE" == "True" ] ;then
	LINEAR_HALLMARK_CONTIGS=$( find . -maxdepth 1 -type f -regextype sed -regex ".*_vs[0-9]\{1,2\}.fna" )
else
	LINEAR_HALLMARK_CONTIGS=$( find . -maxdepth 1 -type f -regextype sed -regex ".*.fna" )
fi


if [ -n "$LINEAR_HALLMARK_CONTIGS" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running BLASTX, annotate linear contigs " $MDYT
	echo "$LINEAR_HALLMARK_CONTIGS" | sed 's/.fna//g' | xargs -n 1 -I {} -P $CPU -t blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -threshold 21 -word_size 5 -num_threads 1 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190925 -query {}.fna -out {}.tax_guide.blastx.out >/dev/null 2>&1
	echo "$LINEAR_HALLMARK_CONTIGS" | while read nucl_fa ; do
		if [ ! -s "${nucl_fa%.fna}.tax_guide.blastx.out" ]; then
			echo "No homologues found" > ${nucl_fa%.fna}.tax_guide.blastx.out ;
		elif grep -i -q "circovir\|genomovir\|geminivir\|nanovir\|redondovir\|bacilladnavir\|smacovir" ${nucl_fa%.fna}.tax_guide.blastx.out ; then 
			EVALUE=$( head -n1 "${nucl_fa%.fna}.tax_guide.blastx.out" | cut -f4 ) ; 
			NEW_TAX=$( head -n1 ${nucl_fa%.fna}.tax_guide.blastx.out | awk -v VALUE="$EVALUE" '{if (VALUE>1e-50) { print $0 ; print "CRESS virus" } else { print $0}}' )
			echo "$NEW_TAX" > ${nucl_fa%.fna}.tax_guide.blastx.out ;
			if grep -q "CRESS virus" ${nucl_fa%.fna}.tax_guide.blastx.out ; then
				echo ${nucl_fa%.fna} "is a CRESS virus"
			else
				ktClassifyBLAST -o ${nucl_fa%.fna}.tax_guide.blastx.tab ${nucl_fa%.fna}.tax_guide.blastx.out >/dev/null 2>&1
				taxid=$( tail -n1 ${nucl_fa%.fna}.tax_guide.blastx.tab | cut -f2 )
				efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element Lineage >> ${nucl_fa%.fna}.tax_guide.blastx.out	
				sleep 0.4s
			fi
		elif grep -q "virophage" ${nucl_fa%.fna}.tax_guide.blastx.out ; then
			echo "Virophage" >> ${nucl_fa%.fna}.tax_guide.blastx.out
		elif grep -q "adinto" ${nucl_fa%.fna}.tax_guide.blastx.out ; then
			echo "Adintovirus" >> ${nucl_fa%.fna}.tax_guide.blastx.out
		elif grep -i -q "polinton" ${nucl_fa%.fna}.tax_guide.blastx.out ; then
			echo "Polinton-like virus" >> ${nucl_fa%.fna}.tax_guide.blastx.out
		else
			ktClassifyBLAST -o ${nucl_fa%.fna}.tax_guide.blastx.tab ${nucl_fa%.fna}.tax_guide.blastx.out >/dev/null 2>&1
			taxid=$( tail -n1 ${nucl_fa%.fna}.tax_guide.blastx.tab | cut -f2 )
			efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -element Lineage >> ${nucl_fa%.fna}.tax_guide.blastx.out
			sleep 0.4s
		fi
		if [ ! -s ${nucl_fa%.fna}.tax_guide.blastx.out ] ; then
			echo "No homologues found" > ${nucl_fa%.fna}.tax_guide.blastx.out
		fi
		if grep -i -q "Caudovir\|Ackermannvir\|Herellevir\|Corticovir\|Levivir\|Tectivir\|crAss-like virus\|CrAssphage\|Cyanophage\|Microvir\microphage\|Siphoviridae\|Myoviridae\|phage\|Podovir\|Halovir\|sphaerolipovir\|pleolipovir\|plasmid\|Inovir\|Ampullavir\|Bicaudavir\|Fusellovir\|Guttavir\|Ligamenvir\|Plasmavir\|Salterprovir\|Cystovir" ${nucl_fa%.fna}.tax_guide.blastx.out ; then
			echo ${nucl_fa} >> LIN_seqs_for_phanotate.txt
		else
			echo ${nucl_fa} >> LIN_seqs_for_prodigal.txt
		fi
	done
	if [ -s LIN_seqs_for_phanotate.txt ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running PHANOTATE, annotate linear contigs " $MDYT
		cat LIN_seqs_for_phanotate.txt | sed 's/.fna//g ; s/\.\///g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/PHANOTATE/phanotate.py -f fasta -o {}.phan.fasta {}.fna
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
				ORIG_CONTIG=$( grep ">" ${PHAN%.phan.fasta}.fna | cut -d " " -f2 ) ; 
				AA_SEQ=$( echo "$LINE" | cut -f2 ) ; 
				if echo $AA_SEQ | grep -q "\*" ; then
					INC3=""
				else
					INC3="3primeInc"
				fi
				FAA=${AA_SEQ:0:1}
				if [ "$FAA" != "M" ] && [ $START_BASE -le 3 ]; then
					INC5="5primeInc"
				else
					INC5=""
				fi

				let COUNTER=COUNTER+1 ; 
				echo ">"${ORF_NAME}"_"${COUNTER} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${PHAN%.phan.fasta}.AA.fasta
		done			
	fi
	if [ -s LIN_seqs_for_prodigal.txt ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running Prodigal, annotate linear contigs " $MDYT
		cat LIN_seqs_for_prodigal.txt | sed 's/.fna//g ; s/\.\///g' | xargs -n 1 -I {} -P $CPU prodigal -a {}.prodigal.fasta -i {}.fna -p meta -q >/dev/null 2>&1
		for PROD in *prodigal.fasta ; do
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
				ORIG_CONTIG=$( grep ">" ${PROD%.prodigal.fasta}.fna | cut -d " " -f2 ) ; 
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
				echo ">"${ORF_NAME} "["$START_BASE" - "$END_BASE"]" ${INC5}${INC3} $ORIG_CONTIG ; echo $AA_SEQ ; 
			done > ${PROD%.prodigal.fasta}.AA.fasta
		done
	fi
	TRANS_AAs=$( find . -maxdepth 1 -type f -name "${run_title}*AA.fasta" | sed 's/\.\///g' )
	if [ -n "$TRANS_AAs" ] ; then
		for ROT in $TRANS_AAs ; do 
			bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' $ROT > ${ROT%.fasta}.sorted.fasta
		done
	fi
fi

#-# 4 hhmscan linear contigs
LIN_SORT_AAs=$( find . -maxdepth 1 -type f -name "${run_title}*AA.sorted.fasta" | sed 's/\.\///g' )
if [ -n "$LIN_SORT_AAs" ] ; then
	cat $( find . -maxdepth 1 -type f -name "${run_title}*AA.sorted.fasta" ) > all_LIN_sort_genome_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_LIN_sort_genome_proteins.AA.fasta | wc -l | bc )
	AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
	if [ $AA_SEQS_PER_FILE = 0 ] ; then
		AA_SEQS_PER_FILE=1
	fi
	awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_LIN_sort_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_LIN_sort_genome_proteins.AA.fasta
	SPLIT_AA_DTR_sort=$( find . -maxdepth 1 -type f -name "SPLIT_LIN_sort_GENOME_AA_*.fasta" )
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running hmmscan1, annotating linear contigs " $MDYT
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
	HMM_REP_NUMEBR=$( find . -maxdepth 1 -type f -name "SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan_replicate.out" | wc -l )
	if [[ $FOR_PLASMIDS = "True" ]]; then
		if [ $HMM_REP_NUMEBR -gt 0 ] ; then
			cat SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan.out SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_LIN_COMBINED.AA.hmmscan.sort.out
		else
			cat SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan.out | grep -v "^#\|plasmid_clust" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_LIN_COMBINED.AA.hmmscan.sort.out
		fi
	else
		if [ $HMM_REP_NUMEBR -gt 0 ] ; then
			cat SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan.out SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan_replicate.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_LIN_COMBINED.AA.hmmscan.sort.out
		else
			cat SPLIT_LIN_sort_GENOME_AA_*AA.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_LIN_COMBINED.AA.hmmscan.sort.out
		fi
	fi		
	if [ -s SPLIT_LIN_COMBINED.AA.hmmscan.sort.out ] ; then
		cut -f3 SPLIT_LIN_COMBINED.AA.hmmscan.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
			grep "${HIT}_" SPLIT_LIN_COMBINED.AA.hmmscan.sort.out > ${HIT}.AA.hmmscan.sort.out
			grep "${HIT}_" SPLIT_LIN_COMBINED.AA.hmmscan.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan1.txt
			HMMSCAN_NUM=$( cat ${HIT}.AA.called_hmmscan1.txt | wc -l | bc )
			TOT_AA_NUM=$( grep -F ">" ${HIT}.AA.sorted.fasta | wc -l | bc )
			if [ $HMMSCAN_NUM -lt $TOT_AA_NUM ] ; then
				grep -v -f ${HIT}.AA.called_hmmscan1.txt ${HIT}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.AA.no_hmmscan1.fasta
			else
				echo " " > ${HIT}.AA.no_hmmscan1.fasta
			fi
		done
	fi
	for LIN in $LIN_SORT_AAs ; do 
		if [ ! -s ${LIN%.AA.sorted.fasta}.AA.no_hmmscan1.fasta ] ; then
			cp $LIN ${LIN%.AA.sorted.fasta}.AA.no_hmmscan1.fasta
		fi
	done
	DTR_AA_FOR_HMM2=$( find . -maxdepth 1 -type f -name "${run_title}*AA.no_hmmscan1.fasta" )
	if [ -n "$DTR_AA_FOR_HMM2" ] ; then
		cat $( find . -maxdepth 1 -type f -name "${run_title}*AA.no_hmmscan1.fasta" ) > all_LIN_HMM2_proteins.AA.fasta
		TOTAL_AA_SEQS=$( grep -F ">" all_LIN_HMM2_proteins.AA.fasta | wc -l | bc )
		if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
			AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
			if [ $AA_SEQS_PER_FILE = 0 ] ; then
				AA_SEQS_PER_FILE=1
			fi
			awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_LIN_HMM2_GENOME_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_LIN_HMM2_proteins.AA.fasta
			SPLIT_LIN_HMM2=$( find . -maxdepth 1 -type f -name "SPLIT_LIN_HMM2_GENOME_AA_*.fasta" )
			MDYT=$( date +"%m-%d-%y---%T" )
			echo "time update: running hmmscan2, annotating linear contigs " $MDYT
			echo "$SPLIT_LIN_HMM2" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --tblout {}.AA.hmmscan2.out --cpu 1 -E 1e-8 --noali ${CENOTE_SCRIPT_DIR}/hmmscan_DBs/useful_hmms_baits_and_not2a {}.fasta >/dev/null 2>&1
			cat SPLIT_LIN_HMM2_GENOME_AA_*AA.hmmscan2.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > SPLIT_LIN_HMM2_COMBINED.AA.hmmscan2.sort.out
		else
			echo "no AA seqs for hmmscan2, annotating linear contigs"
		fi
	fi
	if [ -s SPLIT_LIN_HMM2_COMBINED.AA.hmmscan2.sort.out ] ; then
		cut -f3 SPLIT_LIN_HMM2_COMBINED.AA.hmmscan2.sort.out | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read HIT ; do
			grep "${HIT}_" SPLIT_LIN_HMM2_COMBINED.AA.hmmscan2.sort.out > ${HIT}.AA.hmmscan2.sort.out
			grep "${HIT}_" SPLIT_LIN_HMM2_COMBINED.AA.hmmscan2.sort.out | sort -u -k3,3 | cut -f3 | sed 's/\(.*\)/\1 /' > ${HIT}.AA.called_hmmscan2.txt
			HMMSCAN_NUM=$( cat ${HIT}.AA.called_hmmscan2.txt | wc -l | bc )
			TOT_AA_NUM=$( grep -F ">" ${HIT}.AA.no_hmmscan1.fasta | wc -l | bc )
			if [ $HMMSCAN_NUM -lt $TOT_AA_NUM ] ; then
				grep -v -f ${HIT}.AA.called_hmmscan2.txt ${HIT}.AA.no_hmmscan1.fasta | grep -A1 ">" | sed '/--/d' > ${HIT}.AA.no_hmmscan2.fasta
			else
				echo " " > ${HIT}.AA.no_hmmscan2.fasta
			fi
		done

	fi
	for LIN in $LIN_SORT_AAs ; do 
		if [ ! -s ${LIN%.AA.sorted.fasta}.AA.no_hmmscan2.fasta ] ; then
			cp ${LIN%.AA.sorted.fasta}.AA.no_hmmscan1.fasta ${LIN%.AA.sorted.fasta}.AA.no_hmmscan2.fasta
		fi
	done
	for ROT_AAs in $LIN_SORT_AAs ; do
		echo ">Feature "${ROT_AAs%.AA.sorted.fasta}" Table1" > ${ROT_AAs%.AA.sorted.fasta}.SCAN.tbl
		CALL_ALL_HMM=$( find . -maxdepth 1 -type f -regextype sed -regex "./${ROT_AAs%.AA.sorted.fasta}.*called_hmmscan.*txt" )
		if [ -n "$CALL_ALL_HMM" ] ; then
			cat $( find . -maxdepth 1 -type f -regextype sed -regex "./${ROT_AAs%.AA.sorted.fasta}\..*called_hmmscan.*txt" ) > ${ROT_AAs%.AA.sorted.fasta}.all_called_hmmscans.txt
			if [ -s ${ROT_AAs%.AA.sorted.fasta}.all_called_hmmscans.txt ] ; then
				cat ${ROT_AAs%.AA.sorted.fasta}.all_called_hmmscans.txt | sed 's/ $//g' | while read LINE ; do 
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
					echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""protein motif:$INFERENCEH" >> ${ROT_AAs%.AA.sorted.fasta}.SCAN.tbl ;
				done
			else
				echo "${ROT_AAs%.rotate.AA.sorted.fasta} (linear contig) has no hits in hmmscan1 or hmmscan2 databases"
			fi
		fi
	done
fi


#- blastn all linear seqs
if [ -n "$LINEAR_HALLMARK_CONTIGS" ] && [ $handle_knowns == "blast_knowns" ] ; then
	if [ -s ${BLASTN_DB}.nsq ] || [ -s ${BLASTN_DB}.1.nsq ] || [ -s ${BLASTN_DB}.01.nsq ] ; then
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running BLASTN, linear contigs " $MDYT

		### new blastn
		echo "$LINEAR_HALLMARK_CONTIGS" | sed 's/.fna//g' | xargs -n 1 -I {} -P $CPU blastn -task megablast -query {}.fna -db ${BLASTN_DB} -outfmt '6 std qlen slen' -max_target_seqs 100 -perc_identity 90 -num_threads 1 -word_size 26 -evalue 1e-20 -out {}.blastn.out >/dev/null 2>&1
		for nucl_fa in $LINEAR_HALLMARK_CONTIGS ; do
			if [ -s "${nucl_fa%.fna}.blastn.out" ]; then
				python ${CENOTE_SCRIPT_DIR}/anicalc/anicalc.py -i ${nucl_fa%.fna}.blastn.out -o ${nucl_fa%.fna}.blastn_anicalc.out
				awk '{OFS="\t"}{FS="\t"}{ if (NR==1) {print $1, $2, $4, $5} else if ($4>=95 && $5>=85) {print $1, $2, $4, $5}}' ${nucl_fa%.fna}.blastn_anicalc.out | head -n2 > ${nucl_fa%.fna}.blastn_intraspecific.out
			fi
			if [ -s "${nucl_fa%.fna}.blastn_intraspecific.out" ]; then
				INTRA_LINES=$( cat ${nucl_fa%.fna}.blastn_intraspecific.out | wc -l | bc )	
				if [ "$INTRA_LINES" -ge 2 ] ; then
					ktClassifyBLAST -o ${nucl_fa%.fna}.tax_guide.blastn.tab ${nucl_fa%.fna}.blastn_intraspecific.out >/dev/null 2>&1
					taxid=$( grep -v "qname" ${nucl_fa%.fna}.tax_guide.blastn.tab | tail -n+2 | head -n1 | cut -f2 )
					efetch -db taxonomy -id $taxid -format xml | xtract -pattern Taxon -tab "\n" -element Lineage ScientificName > ${nucl_fa%.fna}.tax_guide.blastn.out
					sleep 1s
					if [ !  -z "${nucl_fa%.fna}.tax_guide.blastn.out" ] ; then

						if grep -i -q "virus\|viridae\|virales\|Circular-genetic-element\|Circular genetic element\|plasmid\|phage" ${nucl_fa%.fna}.tax_guide.blastn.out ; then
							cp ${nucl_fa%.fna}.tax_guide.blastn.out ${nucl_fa%.fna}.tax_guide.KNOWN_VIRUS.out
						else 
							cp ${nucl_fa%.fna}.tax_guide.blastn.out ${nucl_fa%.fna}.tax_guide.CELLULAR.out
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

PROTEIN_NO_HMMSCAN2=$( find . -maxdepth 1 -type f -name "*.AA.no_hmmscan2.fasta" | sed 's/\.\///g' )

if [ -n "$PROTEIN_NO_HMMSCAN2" ]; then

	cat $( find . -maxdepth 1 -type f -name "*.AA.no_hmmscan2.fasta" ) > all_LIN_rps_proteins.AA.fasta
	TOTAL_AA_SEQS=$( grep -F ">" all_LIN_rps_proteins.AA.fasta | wc -l | bc )
	if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_LIN_RPS_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_LIN_rps_proteins.AA.fasta
		SPLIT_LIN_AA_RPS=$( find . -maxdepth 1 -type f -name "SPLIT_LIN_RPS_AA_*.fasta" )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "time update: running RPSBLAST, annotating linear contigs " $MDYT
		echo "$SPLIT_LIN_AA_RPS" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t rpsblast -evalue 1e-4 -num_descriptions 5 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/cdd_rps_db/Cdd -seg yes -query {}.fasta -line_length 200 -out {}.rpsb.out >/dev/null 2>&1
		cat *rpsb.out | awk '{ if ($0 ~ /^>/) {printf $0 ; getline; print $0} else { print $0}}' > COMBINED_RESULTS.rotate.AA.rpsblast.out
		perl ${CENOTE_SCRIPT_DIR}/rpsblastreport2tbl_mt_annotation_pipe_biowulf.pl ;
		grep "protein_id	" COMBINED_RESULTS.NT.tbl | sed 's/.*protein_id	lcl|//g' | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read CONTIG ; do
			#echo "${CONTIG}"
			echo ">Feature ${CONTIG} Table1" > ${CONTIG}.NT.tbl

			grep -A2 -B1 "${CONTIG}_" COMBINED_RESULTS.NT.tbl | sed '/--/d' >> ${CONTIG}.NT.tbl
		done
	else
		echo "no AA seqs for RPSBLAST, annotating linear contigs"
	fi
fi

# Detecting any tRNAs and making a tbl addenum file
if [ -n "$LINEAR_HALLMARK_CONTIGS" ]; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "time update: running tRNAscan-SE " $MDYT
	for GENOME_NAME in $LINEAR_HALLMARK_CONTIGS ; do
		tRNAscan-SE -Q -G -o ${GENOME_NAME}.trnascan-se2.txt ${GENOME_NAME} >/dev/null 2>&1
		
		if grep -q "${GENOME_NAME%.fna}" ${GENOME_NAME}.trnascan-se2.txt ;then
			#echo "$(tput setaf 5) "$GENOME_NAME" was found to encode tRNA(s); making .tbl file $(tput sgr 0)"

			grep "${GENOME_NAME%.fna}" $GENOME_NAME.trnascan-se2.txt | while read LINE ; do 
				TRNA_START=$( echo $LINE | cut -d " " -f3 ) ; 
				TRNA_END=$( echo $LINE | cut -d " " -f4 ) ; 
				TRNA_NUMBER=$( echo $LINE | cut -d " " -f2 ) ; 
				TRNA_TYPE=$( echo $LINE | cut -d " " -f5 ) ; 
				TRNA_SCORE=$( echo $LINE | cut -d " " -f9 ) ; 
				echo -e "$TRNA_START\t""$TRNA_END\t""tRNA\n""\t\t\tgene\t""${GENOME_NAME%.fna}""_tRNA$TRNA_NUMBER\n""\t\t\tproduct\t""tRNA-$TRNA_TYPE\n""\t\t\tinference\t""tRNAscan-SE score:$TRNA_SCORE" >> ${GENOME_NAME%.fna}.trna.tbl; 
			done
		fi
	done
fi


if [ -n "$LINEAR_HALLMARK_CONTIGS" ]; then
	for nucl_fa in $LINEAR_HALLMARK_CONTIGS ; do
		if [ -s ${nucl_fa%.fna}.NT.tbl ] && [ -s ${nucl_fa%.fna}.SCAN.tbl ] && [ -s ${nucl_fa%.fna}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fna}.BLASTP.tbl > ${nucl_fa%.fna}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
			grep -v -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' ${nucl_fa%.fna}.NT.tbl | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' >> ${nucl_fa%.fna}.int.tbl ;
			echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
			cat ${nucl_fa%.fna}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.SCAN.tbl ] && [ -s ${nucl_fa%.fna}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fna}.BLASTP.tbl > ${nucl_fa%.fna}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
			cat ${nucl_fa%.fna}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.NT.tbl ] && [ -s ${nucl_fa%.fna}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fna}.BLASTP.tbl > ${nucl_fa%.fna}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
			grep -v -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' ${nucl_fa%.fna}.NT.tbl | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' >> ${nucl_fa%.fna}.int.tbl ;
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.NT.tbl ] && [ -s ${nucl_fa%.fna}.SCAN.tbl ] ; then
			cat ${nucl_fa%.fna}.NT.tbl > ${nucl_fa%.fna}.int.tbl
			echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
			cat ${nucl_fa%.fna}.SCAN.tbl | grep -v ">Feature" >> ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.NT.tbl ] ; then
			cat ${nucl_fa%.fna}.NT.tbl > ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.SCAN.tbl ] ; then

			cat ${nucl_fa%.fna}.SCAN.tbl >> ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
		elif [ -s ${nucl_fa%.fna}.BLASTP.tbl ] ; then
			cat ${nucl_fa%.fna}.BLASTP.tbl > ${nucl_fa%.fna}.int.tbl
			if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
				echo -e "\n" >> ${nucl_fa%.fna}.int.tbl
				cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
			fi
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
			GENOME_LENGTH=$( bioawk -c fastx '{print length($seq)}' ${feat_tbl3%.int.tbl}.fna )
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
		grep -f ${blastp_tbl1%.int2.tbl}.for_hhpred.txt -A1 ${blastp_tbl1%.int2.tbl}.AA.sorted.fasta | sed '/--/d' > ${blastp_tbl1%.int2.tbl}.blast_hypo.fasta ;
		csplit -z ${blastp_tbl1%.int2.tbl}.blast_hypo.fasta '/>/' '{*}' --prefix=${blastp_tbl1%.int2.tbl} --suffix-format=%02d.for_hhpred.fasta >/dev/null 2>&1 
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
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/hh-suite/build/src/hhsearch -i {}.for_hhpred.fasta -d $PDB_HHSUITE -d $PFAM_HHSUITE -d $CD_HHSUITE -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done		
	elif [[ $HHSUITE_TOOL = "hhblits" ]] ; then
		echo "$dark_orf_list" | sed 's/.for_hhpred.fasta//g' | xargs -n 1 -I {} -P $CPU ${CENOTE_SCRIPT_DIR}/hh-suite/build/src/hhblits -i {}.for_hhpred.fasta -d $PDB_HHSUITE -d $PFAM_HHSUITE -d $CD_HHSUITE -o {}.out.hhr -cpu 1 -maxmem 1 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 >/dev/null 2>&1
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	else
		echo "$(tput setaf 5) Valid option for HHsuite tool (i.e. hhsearch or hhblits) was not provided. Skipping step for "$dark_orf" $(tput sgr 0)"
		for dark_orf in $dark_orf_list ; do	
			cat $dark_orf >> ${dark_orf%*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
			rm -f $dark_orf
		done
	fi
	cat_list=$( find . -maxdepth 1 -type f -name "*.out.hhr" )
	if [ -n "$cat_list" ] ; then
		cat *out.hhr > ${run_title}.rotate.out_all.hhr
		rm -f *out.hhr
	fi
fi

rm -f *[0-9].AA.fasta

perl ${CENOTE_SCRIPT_DIR}/hhpredreport2tbl_mt_annotation_pipe_biowulf1_gjs_edits.pl
if [ -s ${run_title}.HH.tbl ] ; then
	grep "protein_id	" ${run_title}.HH.tbl | sed 's/.*protein_id	lcl|//g' | sed 's/[^_]*$//' | sed 's/\(.*\)_/\1/' | sort -u | while read CONTIG ; do
		echo ">Feature ${CONTIG} Table1" > ${CONTIG}.HH.tbl
		grep -A2 -B1 "${CONTIG}_" ${run_title}.HH.tbl | sed '/--/d' >> ${CONTIG}.HH.tbl
	done
	mv ${run_title}.HH.tbl ${run_title}.combined_hhsuite.tbl
fi

HH_TBL=$( find . -maxdepth 1 -type f -name "*.HH.tbl" )
if [ -n "$HH_TBL" ] ; then
	for HH_tbl1 in $HH_TBL ; do 
		sed 's/OS=.*//g; s/ ;//g; s/similar to AA sequence:UniProtKB:>\([0-9][A-Z].*\)/protein motif:PDB:\1/g; s/UniProtKB:>tr|.*|\(.\)/UniProtKB:\1/g; s/similar to AA sequence:UniProtKB:>\([a-z].*\)/protein motif:Scop:\1/g; s/similar to AA sequence:UniProtKB:>\(PF.*\)/protein motif:PFAM:\1/g; s/ is .*//g; s/ are .*//g' $HH_tbl1 | sed '/product/ s/; [a-zA-Z0-9_]\{1,20\}//g; s/;.*//g' > ${HH_tbl1%.HH.tbl}.HH2.tbl
	done
fi

# Combining tbl files from all search results AND fix overlapping ORF module

INT2_TBL=$( find . -maxdepth 1 -type f -name "*.int2.tbl" )

if [ -n "$INT2_TBL" ] ; then
	echo "$(tput setaf 5) Combining tbl files from all search results AND fix overlapping ORF module, linear contigs $(tput sgr 0)"
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
		GENOME_LENGTH=$( bioawk -c fastx '{print length($seq)}' ${feat_tbl4%.int2.tbl}.fna )
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
else
	echo "int2.tbl not found"
fi
### insert comb3.tbl edits
COMB3_TBL=$( find . -maxdepth 1 -type f -name "*.comb3.tbl" )
if [ -n "$COMB3_TBL" ] ; then
	for comb3 in $COMB3_TBL ; do
		## tRNA overlap
		if grep -q "[0-9]	tRNA" $comb3 ; then
			grep "[0-9]	tRNA" $comb3 | awk '{OFS="\t"}{FS="\t"}{ if ($1<$2) {print "name", $1, $2, "fwd"} else {print "name", $2, $1, "rev"}}' > ${comb3%.comb3.tbl}.tRNA.bed
		fi
		if grep -q "	CDS" $comb3 ; then
			grep "	CDS" $comb3 | awk '{OFS="\t"}{FS="\t"}{ if ($1<$2) {print "name", $1, $2, "fwd"} else {print "name", $2, $1, "rev"}}' > ${comb3%.comb3.tbl}.CDS.bed
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
		if grep -q "5primeInc\|3primeInc" ${comb3%.comb3.tbl}.AA.sorted.fasta ; then
			grep "5primeInc\|3primeInc" ${comb3%.comb3.tbl}.AA.sorted.fasta | while read INCOMPLETE ; do
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
	echo "finalizing taxonomy for linear contigs"
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

				grep -A1 "$TAX_ORF " ${feat_tbl2%.comb3.tbl}.AA.sorted.fasta | sed '/--/d' > ${feat_tbl2%.comb3.tbl}.tax_orf.fasta
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

MDYT=$( date +"%m-%d-%y---%T" )
echo "time update: finished annotating linear contigs " $MDYT