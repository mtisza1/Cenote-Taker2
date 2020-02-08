#!/bin/bash

echo "$(tput setaf 5) Starting Sketch of putatively non-viral contigs $(tput sgr 0)"

cd $base_directory/$run_title
mkdir noncircular_non_viral_domain_contig_sketches
cd noncircular_contigs/


if [ -z noncircular_non_viral_domains_contigs.fna ] ; then
	echo "$(tput setaf 5)There are no contigs in the noncircular_non_viral_domains_contigs.fna file, exiting module $(tput sgr 0)"
else
	grep "^>" noncircular_non_viral_domains_contigs.fna | sed 's/>//g' | while read LINE ; do
		CONTIG_NAME=$( echo $LINE | cut -f1 -d " " ) 
		grep -A1 "$LINE" noncircular_non_viral_domains_contigs.fna | sed '/--/d' > ../noncircular_non_viral_domain_contig_sketches/$CONTIG_NAME.fasta ; 
	done


	cd ../noncircular_non_viral_domain_contig_sketches

	NUCL_FILES_LIST=$( ls *.fasta )
	for NUCL_FILES in $NUCL_FILES_LIST ; do
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -find 1 -minsize 240 -sequence $NUCL_FILES -outseq ${NUCL_FILES%.fasta}.AA.fasta ;
	done


	time ls *.AA.fasta | xargs -n 1 -I {} -P $CPU -t rpsblast -evalue 1e-4 -num_descriptions 5 -num_alignments 1 -db ${CENOTE_SCRIPT_DIR}/cdd_rps_db/Cdd -query {} -line_length 100 -out {}.rpsblast1.out

		echo "$(tput setaf 5)RPS-BLAST of putative non-viral contigs complete.$(tput sgr 0)"
		echo " "

	perl ${CENOTE_SCRIPT_DIR}/rpsblastreport2tbl_mt_sketch_contigs1.pl ;

	echo "$(tput setaf 5) Looking for tRNAs in contigs; non-circular/non-ITR contigs without viral domains  $(tput sgr 0)"
	for GENOME_NAME in $NUCL_FILES_LIST ; do
		tRNAscan-SE -Q -G -o $GENOME_NAME.trnascan-se2.txt $GENOME_NAME
		
		if grep -q "${GENOME_NAME%.fasta}" $GENOME_NAME.trnascan-se2.txt ;then
			echo "$(tput setaf 5) "$GENOME_NAME" was found to encode tRNA(s); making .tbl file $(tput sgr 0)"

			grep "${GENOME_NAME%.fasta}" $GENOME_NAME.trnascan-se2.txt | while read LINE ; do 
				TRNA_START=$( echo $LINE | cut -d " " -f3 ) ; 
				TRNA_END=$( echo $LINE | cut -d " " -f4 ) ; 
				TRNA_NUMBER=$( echo $LINE | cut -d " " -f2 ) ; 
				TRNA_TYPE=$( echo $LINE | cut -d " " -f5 ) ; 
				TRNA_SCORE=$( echo $LINE | cut -d " " -f9 ) ; 
				echo -e "$TRNA_START\t""$TRNA_END\t""tRNA\n""\t\t\tgene\t""$GENOME_NAME""_tRNA$TRNA_NUMBER\n""\t\t\tproduct\t""tRNA-$TRNA_TYPE\n""\t\t\tinference\t""tRNAscan-SE score:$TRNA_SCORE" >> ${GENOME_NAME%.fasta}.trna.tbl; 
			done
		fi
	done

	for GENOME_NAME in $NUCL_FILES_LIST ; do
		if [ -s ${GENOME_NAME%.fasta}.NT.tbl ] ; then
			cat ${GENOME_NAME%.fasta}.NT.tbl > ${GENOME_NAME%.fasta}.final.tbl
			if [ -s ${GENOME_NAME%.fasta}.trna.tbl ] ; then
				echo -e "\n" >> ${GENOME_NAME%.fasta}.final.tbl
				cat ${GENOME_NAME%.fasta}.trna.tbl >> ${GENOME_NAME%.fasta}.final.tbl
			fi
		else
			echo "$(tput setaf 5) "$GENOME_NAME" did not have an associated .tbl file $(tput sgr 0)"
		fi
	done

	echo "$(tput setaf 5) Making .gff files for each annotated sequence $(tput sgr 0)"

	for feat_tbl2 in *.final.tbl ; do
		if [ -s ${feat_tbl2%.final.tbl}.gtf ] ; then
			rm ${feat_tbl2%.final.tbl}.gtf
		fi
		grep "^[0-9]" -A3 $feat_tbl2 | sed '/--/d' | sed 's/ /_/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n		//g' | while read LINE ; do
			if echo $LINE | grep -q "CDS" ; then
				GENOME=${feat_tbl2%.final.tbl}
				FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
				FEAT_START=$( echo $LINE | cut -d " " -f1 )
				FEAT_END=$( echo $LINE | cut -d " " -f2 )
				FEAT_NAME=$( echo $LINE | cut -d " " -f7 )
				FEAT_ATT=$( echo $LINE | cut -d " " -f9 )
				FEAT_ID=$( echo $LINE | cut -d " " -f5 )
			elif echo $LINE | grep -q "repeat_region" ; then
				GENOME=${feat_tbl2%.final.tbl}
				FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
				FEAT_START=$( echo $LINE | cut -d " " -f1 )
				FEAT_END=$( echo $LINE | cut -d " " -f2 )
				FEAT_NAME="ITR"
				FEAT_ATT="ITR"
				FEAT_ID="ITR"	
			fi


			echo -e "$GENOME\t""Cenote-Taker\t""$FEAT_TYPE\t""$FEAT_START\t""$FEAT_END\t"".\t"".\t"".\t""gene_id \"$FEAT_ID\"; gene_name \"$FEAT_NAME\"; gene_inference \"$FEAT_ATT\"" >> ${feat_tbl2%.final.tbl}.gtf
		done
	done

	# Making directory for sequin generation
	if [ ! -d "sequin_directory" ]; then
		mkdir sequin_directory
	fi

	for feat_tbl2 in *.final.tbl ; do 
		file_core=${feat_tbl2%.final.tbl}
		file_numbers=$( echo ${file_core: -3} | sed 's/[a-z]//g' | sed 's/[A-Z]//g' )
		vir_name="Undetermined"
		rand_id=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )
		seq_name1=$( head -n1 ${feat_tbl2%.final.tbl}.fasta | sed 's/>//g; s/|.*//g' | cut -d " " -f2 )

		if [ -s $feat_tbl2 ]; then 
			cp $feat_tbl2 sequin_directory/${feat_tbl2%.final.tbl}.tbl ; 
			bioawk -v contig_name="$seq_name1" -v srr_var="$srr_number" -v headername="$vir_name" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= "contig_name" ; taxonomy was not determined ] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.final.tbl}.fasta > sequin_directory/${feat_tbl2%.final.tbl}.fsa ;	
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
		fi
	done


	if [[ $DATA_SOURCE = "tpa_assembly" ]] ;then
		${CENOTE_SCRIPT_DIR}/linux64.tbl2asn -V vb -j "[keyword=TPA:assembly]" -t $base_directory/$template_file -X C -p sequin_directory/ ;
	else
		${CENOTE_SCRIPT_DIR}/linux64.tbl2asn -V vb -t $base_directory/$template_file -X C -p sequin_directory/ ;
	fi
fi

cd ..

echo "$(tput setaf 3) FINISHED SKETCHING CONTIGS WITHOUT ANY VIRAL DOMAINS OR CIRCULARITY OR ITRS USING RPS-BLAST $(tput sgr 0)"



