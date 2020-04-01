#!/bin/bash

# This is the Cenote-Taker script for quick annotation of known circular sequences

cd circles_of_known_viruses
	

known_fastas=$( ls *.fna )
for nucl_fa in $known_fastas ; do
	echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"
	ACCESSIONQ=$( cut -d "|" -f4 ${nucl_fa%.known_species.fna}.blastn.out ) ; 
	efetch -db nuccore -id $ACCESSIONQ -format fasta | bioawk -c fastx '{ print ">"$name ; print substr($seq,1,300) }' > ${nucl_fa%.fna}.starting_orf.1.fa
	if [ -s ${nucl_fa%.fna}.starting_orf.1.fa ] ; then
		circlator fixstart --genes_fa ${nucl_fa%.fna}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fna}.rotate ;
	else
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
			cp $nucl_fa ${nucl_fa%.fna}.no_100AA_ORFs.fasta
		fi
	fi

	taxid=$( tail -n1 ${nucl_fa%.known_species.fna}.tax_guide.blastn.tab | cut -f2 )
	if grep -q "Inovir" ${nucl_fa%.known_species.fna}.tax_guide.blastn.out ; then
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 90 -sequence ${nucl_fa%.fna}.rotate.fasta -outseq ${nucl_fa%.fna}.rotate.AA.fasta ;
		echo "getorf for AA sequence done"
	else
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 150 -sequence ${nucl_fa%.fna}.rotate.fasta -outseq ${nucl_fa%.fna}.rotate.AA.fasta ;
		echo "getorf for AA sequence done"

	fi
	bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' ${nucl_fa%.fna}.rotate.AA.fasta > ${nucl_fa%.fna}.rotate.AA.sorted.fasta
	echo "starting BLASTP of "${nucl_fa%.fna}.rotate.AA.sorted.fasta
	blastp -evalue 1e-10 -num_descriptions 5 -num_threads $CPU -num_alignments 1 -db /fdb/blastdb/nr -query ${nucl_fa%.fna}.rotate.AA.sorted.fasta -out ${nucl_fa%.known_species.fna}.rotate.blastp.out ;

done 
perl ${CENOTE_SCRIPT_DIR}/blastpreport2tbl_mt_annotation_pipe_biowulf2.pl ;

for GENOME_NAME in $novel_fastas ; do
	if grep -i -q "caudovir\|podovir\|siphovir\|myovir\|ackermannvir\|mimivir\|pandoravir\|pithovir\|marseillevir\|plasmid\|" ${GENOME_NAME%.known_species.fna}.tax_guide.blastn.out ; then
		tRNAscan-SE -Q -G -o $GENOME_NAME.trnascan-se2.txt ${GENOME_NAME}
	fi
	
	if grep -q "$GENOME_NAME" $GENOME_NAME.trnascan-se2.txt ;then

		grep "$GENOME_NAME" $GENOME_NAME.trnascan-se2.txt | while read LINE ; do 
			TRNA_START=$( echo $LINE | cut -d " " -f3 ) ; 
			TRNA_END=$( echo $LINE | cut -d " " -f4 ) ; 
			TRNA_NUMBER=$( echo $LINE | cut -d " " -f2 ) ; 
			TRNA_TYPE=$( echo $LINE | cut -d " " -f5 ) ; 
			TRNA_SCORE=$( echo $LINE | cut -d " " -f9 ) ; 
			echo "$TRNA_START\t""$TRNA_END\t""tRNA\n""\t\t\ttrna_id\t""lcl|""$GENOME_NAME""_$TRNA_NUMBER\n""\t\t\tproduct\t""tRNA-$TRNA_TYPE\n""\t\t\tinference\t""tRNAscan-SE score:$TRNA_SCORE" >> ${GENOME_NAME%.fna}.trna.tbl; 
		done
	fi
done

for nucl_fa in $known_fastas ; do
	cat ${nucl_fa%.known_species.fna}.BLASTP.tbl > ${nucl_fa%.known_species.fna}.int.tbl
	if [ -s ${nucl_fa%.fna}.trna.tbl ] ; then
		cat ${nucl_fa%.fna}.trna.tbl >> ${nucl_fa%.fna}.int.tbl
	fi
done
for feat_tbl3 in *.int.tbl ; do
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
	sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g' $feat_tbl3 | sed '/--/d' > ${feat_tbl3%.int.tbl}.comb2.tbl ; 
	grep -e 'hypothetical protein' -B2 ${feat_tbl3%.int.tbl}.comb2.tbl | grep "^[0-9]" | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.int.tbl}.hypo_start_stop.txt ;

	cat "${feat_tbl3%.int.tbl}.hypo_start_stop.txt" | while read liney ; do
		loc_start=$( echo $liney | cut -d " " -f1 )
		loc_end=$( echo $liney | cut -d " " -f2 )
		loc1_start=$( echo " " "$loc_start" " ")
		if grep -q "$loc1_start" ${feat_tbl3%.int.tbl}.used_positions.txt ; then 
			echo "$loc1_start"
			if [[ "$loc_end" -gt "$loc_start" ]]; then
				f_end=$(( $loc_end + 1 ))
				f1_end=$( echo " " "$f_end" " ")
				echo "$f1_end"
				if grep -q "$f1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then
					echo "$loc_end" "end"
					echo "$liney" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
				fi
			else
				r_end=$(( $loc_end - 1 ))
				r1_end=$( echo " " "$r_end" " ")

				echo "$r1_end"
				if grep -q "$r1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then
					echo "$loc_end" "end"
					echo "$liney" >> ${feat_tbl3%.int.tbl}.remove_hypo.txt
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

# Making directory for sequin generation
if [ ! -d "sequin_directory" ]; then
	mkdir sequin_directory
fi

for feat_tbl2 in *.int2.tbl ; do 
	file_core=${feat_tbl2%.int2.tbl}
	echo $file_core
	file_numbers=$( echo ${file_core: -3} | sed 's/[a-z]//g' | sed 's/[A-Z]//g' )
	echo $file_numbers
	tax_guess=$( tail -n1 ${feat_tbl2%.int2.tbl}.tax_guide.blastn.out ) ; 

	vir_name=$( cat ${feat_tbl2%.int2.tbl}.tax_guide.blastn.out | tail -n1 )
	rand_id=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )
	fsa_head=$( echo $vir_name " strain ct"${rand_id}${file_numbers} )
	perc_id=$( head -n1 ${feat_tbl2%.int2.tbl}.tax_guide.blastn.out | sed 's/ /-/g' | awk '{FS="\t"; OFS="\t"} {print $2" "$3}' | sed 's/-/ /g' ) ;

	cp $feat_tbl2 sequin_directory/${feat_tbl2%.int2.tbl}.tbl ; 
	bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.int2.tbl}.known_species.rotate.fasta > sequin_directory/${feat_tbl2%.int2.tbl}.fsa ;	
done

for nucl_fa in $known_fastas ; do
	if echo $ASSEMBLER | grep -i -q "spades" ; then
		coverage=$( head -n1 $nucl_fa | cut -d " " -f 2 | cut -d "_" -f 6 | sed 's/|.*//g' ) 
		echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Assembly Method	" $ASSEMBLER >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Genome Coverage	"$coverage"x" >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Sequencing Technology	Illumina" >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Annotation Method	Cenote-Taker2" >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
	else
		echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Assembly Method	" $ASSEMBLER >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Sequencing Technology	Illumina" >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;
		echo "Annotation Method	Cenote-Taker2" >> sequin_directory/${nucl_fa%.known_species.fna}.cmt ;

	fi
done

${CENOTE_SCRIPT_DIR}/linux64.tbl2asn -V vb -t $base_directory/$template_file -X C -p sequin_directory/ ;

for fsa_file in sequin_directory/*.fsa ; do
	fsa_name2=$( echo ${fsa_file#sequin_directory/} ) ; 
	fsa_name3=$( echo ${fsa_name2%.fsa} | sed 's/.PLASMID//g' )
	seq_name1=$( head -n1 $fsa_name3.known_species.fna | sed 's/>//g; s/|.*//g' | cut -d " " -f2 )
	sed " 1 s/note= closest relative/note= $seq_name1 ; closest relative/" $fsa_file > $fsa_file.temp
	mv $fsa_file.temp $fsa_file
done

${CENOTE_SCRIPT_DIR}/linux64.tbl2asn -V vb -t $base_directory/$template_file -X C -p sequin_directory/ ;

rm *.all_start_stop.txt *.bad_starts.txt *.comb.tbl *.comb2.tbl *.good_start_orfs.txt *.hypo_start_stop.txt *.nucl_orfs.fa *.remove_hypo.txt *.log *.promer.contigs_with_ends.fa *.promer.promer *.out.hhr *.starting_orf.1.fa *.starting_orf.txt *.used_positions.txt
echo "$(tput setaf 3) DONE WITH KNOWN VIRUSES $(tput sgr 0)"

cd ..
