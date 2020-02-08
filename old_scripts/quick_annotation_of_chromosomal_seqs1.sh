#!/bin/bash

# This is the Cenote-Taker script for analysis and annotation of circular sequences highly similar to cellular sequences

cd circles_of_chromosomal_elements

known_fastas=$( ls *.fna )
for nucl_fa in $known_fastas ; do
	echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"
	/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -minsize 240 -find 3 -sequence $nucl_fa -outseq ${nucl_fa%.fna}.nucl_orfs.fa ; 

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
		cp $nucl_fa ${nucl_fa%.fna}.no_100AA_ORFs.fna
	fi
	ktClassifyBLAST -o ${nucl_fa%.fna}.tax_guide.blastn.tab ${nucl_fa%.fna}.blastn.out
	taxid=$( tail -n1 ${nucl_fa%.fna}.tax_guide.blastn.tab | cut -f2 )
	efetch -db taxonomy -id $taxid -format xml | /data/tiszamj/mike_tisza/xtract.Linux -pattern Taxon -element Lineage >> ${nucl_fa%.fna}.tax_guide.blastn.out
	if grep -q "Inovir" ${nucl_fa%.fna}.tax_guide.blastn.out ; then
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 90 -sequence ${nucl_fa%.fna}.rotate.fna -outseq ${nucl_fa%.fna}.rotate.AA.fasta ;
	else
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 240 -sequence ${nucl_fa%.fna}.rotate.fna -outseq ${nucl_fa%.fna}.rotate.AA.fasta ;
	fi
	bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' ${nucl_fa%.fna}.rotate.AA.fasta > ${nucl_fa%.chromosomal.fna}.rotate.AA.sorted.fasta
	hmmscan --tblout ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.out --cpu 56 -E 1e-8 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/baits_edited_1/virus_specific_baits_edited2 ${nucl_fa%.chromosomal.fna}.rotate.AA.sorted.fasta
	echo "$(tput setaf 5)HMMSCAN search of "${nucl_fa%.chromosomal.fna}.rotate.AA.sorted.fasta" complete.$(tput sgr 0)"
	echo " "
	grep -v "^#" ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.out | sed 's/ \+/	/g' | sort -u -k3,3 > ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.sort.out
	cut -f3 ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.sort.out > ${nucl_fa%.chromosomal.fna}.rotate.AA.called_hmmscan.txt ; 
	if [ -s ${nucl_fa%.chromosomal.fna}.rotate.AA.called_hmmscan.txt ] ; then

		grep -v -f ${nucl_fa%.fna}.rotate.AA.called_hmmscan.txt ${nucl_fa%.chromosomal.fna}.rotate.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ../${nucl_fa%.chromosomal.fna}.rotate.no_hmmscan.fasta
		echo ">Feature "${nucl_fa%.chromosomal.fna}" Table1" > ../${nucl_fa%.chromosomal.fna}.SCAN.tbl

		cat ${nucl_fa%.chromosomal.fna}.rotate.AA.called_hmmscan.txt | while read LINE ; do 
			PROTEIN_INFO=$( grep "$LINE" ${nucl_fa%.fna}.rotate.AA.sorted.fasta ) ; 
			START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
			END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
			HMM_INFO=$( grep "$LINE" ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.out | head -n1 | cut -d " " -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g') ; 
			INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
			PROTEIN_NAME=$( echo $HMM_INFO | cut -d " " -f2- ) ; 
			echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""similar to AA sequence:$INFERENCEH" >> ../${nucl_fa%.chromosomal.fna}.SCAN.tbl ; 
		done


		echo "$(tput setaf 5)Guessing taxonomy for sequence "${nucl_fa%.fna}.rotate.fasta" by BLASTX against virus and plasmid protein database.$(tput sgr 0)"
		blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -num_threads 56 -num_alignments 1 -db /data/tiszamj/mike_tisza/auto_annotation_pipeline/blast_DBs/virus_and_plasmid_proteins -query ${nucl_fa%.fna}.rotate.fasta -out ${nucl_fa%.fna}.tax_guide.blastx.out ;
		if [ ! -s "${nucl_fa%.fna}.tax_guide.blastx.out" ]; then
			echo "No homologues found" > ${nucl_fa%.fna}.tax_guide.blastx.out ;
		else
			echo "$(tput setaf 5)"$nucl_fa" likely represents a novel virus or plasmid. Getting hierarchical taxonomy info.$(tput sgr 0)"
			ktClassifyBLAST -o ${nucl_fa%.fna}.tax_guide.blastx.tab ${nucl_fa%.fna}.tax_guide.blastx.out
			taxid=$( tail -n1 ${nucl_fa%.fna}.tax_guide.blastx.tab | cut -f2 )
			efetch -db taxonomy -id $taxid -format xml | /data/tiszamj/mike_tisza/xtract.Linux -pattern Taxon -element Lineage >> ${nucl_fa%.fna}.tax_guide.blastx.out
		fi
		mv ${nucl_fa%.fna}.rotate.fasta ${nucl_fa%.chromosomal.fna}.provirus.fna
		cp ${nucl_fa%.chromosomal.fna}.provirus.fna ../${nucl_fa%.chromosomal.fna}.rotate.fasta
		cp ${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.out ../${nucl_fa%.chromosomal.fna}.rotate.AA.hmmscan.out
		cp ${nucl_fa%.fna}.tax_guide.blastx.out ../${nucl_fa%.fna}.tax_guide.blastx.out
		cp ${nucl_fa%.fna}.tax_guide.blastx.tab ../${nucl_fa%.fna}.tax_guide.blastx.tab


	fi	

done 

cd ..