#!/bin/bash

# This is Cenote-Taker2 for Biowulf.
# See README for inputs. 
# See schematic for how Cenote-Taker2 processes data.


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

# Loading all the modules and environments for biowulf
source /data/tiszamj/conda/etc/profile.d/conda.sh
conda activate cenote_taker1
module load samtools || echo "$(tput setaf 4)unable to load samtools module $(tput sgr 0)"
module load mummer || echo "$(tput setaf 4)unable to load mummer module $(tput sgr 0)"
module load prodigal || echo "$(tput setaf 4)unable to load prodigal module $(tput sgr 0)" 
module load bioawk || echo "$(tput setaf 4)unable to load bioawk module $(tput sgr 0)"
module load blast || echo "$(tput setaf 4)unable to load blast module $(tput sgr 0)" 
module load hhsuite || echo "$(tput setaf 4)unable to load hhsuite module $(tput sgr 0)"
module load bwa || echo "$(tput setaf 4)unable to load BWA module $(tput sgr 0)"
module load edirect || echo "$(tput setaf 4)unable to load edirect module $(tput sgr 0)"
module load kronatools || echo "$(tput setaf 4)unable to load kronatools module $(tput sgr 0)"
module load hmmer || echo "$(tput setaf 4)unable to load hmmer module $(tput sgr 0)"
module load trnascan-se || echo "$(tput setaf 4)unable to load trnascan-se module $(tput sgr 0)"
module load bowtie || echo "$(tput setaf 4)unable to load bowtie2 module $(tput sgr 0)"
module load bbtools || echo "$(tput setaf 4)unable to load bbtools module $(tput sgr 0)"
module load minialign || echo "$(tput setaf 4)unable to load minialign module $(tput sgr 0)"

# Setting input parameters
original_spades_contigs=$1
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
handle_nonviral=${13}
length_cutoff=${14}
virus_domain_db=${15}
MINIMUM_DOMAINS=${16}
handle_knowns=${17}
ASSEMBLER=${18}
MOLECULE_TYPE=${19}
HHSUITE_TOOL=${20}
DATA_SOURCE=${21}
base_directory=$PWD
# Making output folder
if [ ! -d "$run_title" ]; then
	mkdir "$run_title"
fi 

echo "Option -default turned on" 
# Removing contigs under 1000bp and detecting circular contigs
if [ ${original_spades_contigs: -6} == ".fasta" ]; then
	echo "$(tput setaf 5)File with .fasta extension detected, attempting to keep contigs over $length_cutoff nt and find circular sequences with apc.pl$(tput sgr 0)"
	bioawk -v run_var="$run_title" -v contig_cutoff="$length_cutoff" -c fastx '{ if(length($seq) > contig_cutoff) { print ">"run_var NR" "$name; print $seq }}' $original_spades_contigs > ${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta ;
	cd $run_title
	perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/apc_biowulf1.pl -b $run_title ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta ;
	rm apc_aln*
	for fa1 in $run_title*.fa ; do 
		mv $fa1 $run_title${fa1#$run_title.}sta ; 
	done 
elif [ ${original_spades_contigs: -6} == ".fastg" ]; then
	bioawk -v contig_cutoff="$length_cutoff" -c fastx '{ if(length($seq) > contig_cutoff) {print }}' $original_spades_contigs | grep "[a-zA-Z0-9]:\|[a-zA-Z0-9];" | grep -v "':" | awk '{ print ">"$1 ; print $2 }' | sed 's/:.*//g; s/;.*//g' | bioawk -v run_var="$run_title" -c fastx '{ print ">"run_var NR" "$name; print $seq }' > ${original_spades_contigs%.fastg}.over_${length_cutoff}nt.fasta
	cd $run_title
	perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/apc_biowulf1.pl -b $run_title ../${original_spades_contigs%.fastg}.over_${length_cutoff}nt.fasta ;
	rm apc_aln*
	for fa1 in $run_title*.fa ; do 
		mv $fa1 $run_title${fa1#$run_title.}sta ; 
	done 		
else
	echo "$(tput setaf 4)File with .fasta of .fastg extension not detected as first input. Exiting.$(tput sgr 0)" ;
	exit
fi

# Changing to output directory
cd $run_title

# Aligning reads to contigs
echo "$(tput setaf 4)Aligning provided reads to contigs over cutoff to determine coverage. $(tput sgr 0)" 

date +"%m-%d-%y---%T" ;
if [ $F_READS == "-no_reads" ] ; then
	continue
else
	if [ $R_READS == "-unpaired" ] ;then
		mkdir bt2_indices ; 
		echo "$(tput setaf 4)Building Bowtie2 Index. $(tput sgr 0)" 
		bowtie2-build ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta bt2_indices/${run_title}_bt2_index
		echo "$(tput setaf 4)Aligning reads to BowTie2 index. $(tput sgr 0)" 
		bowtie2 -q -p56 -x bt2_indices/${run_title}_bt2_index -U $F_READS -S reads_to_all_contigs_over${length_cutoff}nt.sam --very-fast
	elif [ $R_READS == "-nanopore" ] || [ $R_READS == "-pacbio" ] ; then
		echo "$(tput setaf 4)Aligning Nanopore/Pacbio reads to refernce contigs. $(tput sgr 0)" 
		minialign -t 56 ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta $F_READS > reads_to_all_contigs_over${length_cutoff}nt.sam
	else
		mkdir bt2_indices ; 
		echo "$(tput setaf 4)Building Bowtie2 Index. $(tput sgr 0)" 
		bowtie2-build ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta bt2_indices/${run_title}_bt2_index
		echo "$(tput setaf 4)Aligning reads to BowTie2 index. $(tput sgr 0)" 
		bowtie2 -q -p56 -x bt2_indices/${run_title}_bt2_index -1 $F_READS -2 $R_READS -S reads_to_all_contigs_over${length_cutoff}nt.sam --very-fast
	fi
fi
echo "$(tput setaf 4)Calculating read coverage of each contig with BBTools Pileup. $(tput sgr 0)" 
bbtools pileup in=reads_to_all_contigs_over${length_cutoff}nt.sam out=reads_to_all_contigs_over${length_cutoff}nt.coverage.txt
rm reads_to_all_contigs_over${length_cutoff}nt.sam

# Detecting whether any circular contigs were present
original_fastas=$( ls *.fasta )
# "$(tput setaf 5)$var1$(tput sgr 0)"
if [ -z "$original_fastas" ] ; then
	echo "$(tput setaf 4)No circular fasta files detected. $(tput sgr 0)" 
	#exit
	mkdir noncircular_contigs	
	if [ ${original_spades_contigs: -6} == ".fasta" ]; then
		grep -A1 "^>" ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta | sed 's/--//g' > noncircular_contigs/all_non_circular.fasta
	elif [ ${original_spades_contigs: -6} == ".fastg" ]; then
		grep -A1 "^>" ../${original_spades_contigs%.fastg}.over_${length_cutoff}nt.fasta | sed 's/--//g' > noncircular_contigs/all_non_circular.fasta
	fi			
	if [ -s noncircular_contigs/all_non_circular.fasta ] ; then
		grep "^>" noncircular_contigs/all_non_circular.fasta | sed 's/>//g' | cut -d " " -f1 | while read LINE ; do 
			grep -A1 "$LINE [a-zA-Z]" noncircular_contigs/all_non_circular.fasta > noncircular_contigs/$LINE.fasta ; 
		done
	fi
else
	echo "$(tput setaf 5)Circular fasta file(s) detected$(tput sgr 0)"
	echo " "
# Putting non-circular contigs in a separate directory
	echo "$(tput setaf 4)Putting non-circular contigs in a separate directory $(tput sgr 0)" 

	mkdir noncircular_contigs

	for CIRCLE in $original_fastas ; do
		grep "^>" $CIRCLE | sed 's/|.*//g' >> circular_contigs_spades_names.txt
	done
	if [ ${original_spades_contigs: -6} == ".fasta" ]; then
		grep -v -f circular_contigs_spades_names.txt ../${original_spades_contigs%.fasta}.over_${length_cutoff}nt.fasta | grep -A1 "^>" | sed 's/--//g' > noncircular_contigs/all_non_circular.fasta
	elif [ ${original_spades_contigs: -6} == ".fastg" ]; then
		grep -v -f circular_contigs_spades_names.txt ../${original_spades_contigs%.fastg}.over_${length_cutoff}nt.fasta | grep -A1 "^>" | sed 's/--//g' > noncircular_contigs/all_non_circular.fasta
	fi
	if [ -s noncircular_contigs/all_non_circular.fasta ] ; then
		grep "^>" noncircular_contigs/all_non_circular.fasta | sed 's/>//g' | cut -d " " -f1 | while read LINE ; do 
			grep -A1 "$LINE [a-zA-Z]" noncircular_contigs/all_non_circular.fasta > noncircular_contigs/$LINE.fasta ; 
		done
	fi

fi
# Looking for ITRs in non-circular contigs
echo "$(tput setaf 4)Looking for ITRs in non-circular contigs $(tput sgr 0)" 

cd noncircular_contigs
CONTIGS_NON_CIRCULAR=$( ls *{0..9}.fasta )
echo "$(tput setaf 4)CONTIGS_NON_CIRCULAR variable $(tput sgr 0)"
echo $CONTIGS_NON_CIRCULAR
if [ ! -z "$CONTIGS_NON_CIRCULAR" ] ;then
	for NONCIR in $CONTIGS_NON_CIRCULAR ; do
		LEN_CHECKQ=$( cat $NONCIR | bioawk -c fastx '{ if(length($seq) > 4000) { print $name }}' ) ; 
		if [ ! -z "$LEN_CHECKQ" ] ; then
			/data/tiszamj/mike_tisza/irf307.linux.exe $NONCIR 2 3 5 80 10 40 500000 10000 -d -h 
		fi
	done
	mkdir ../ITR_containing_contigs
	for DAT in *.dat ; do 

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
					mv ${DAT%.2.3.5.80.10.40.500000.10000.dat} ../ITR_containing_contigs/${DAT%.2.3.5.80.10.40.500000.10000.dat}
					### MAKE ITR TBL FILE
					echo "$(tput setaf 4) Making ITR .tbl file $(tput sgr 0)"
					L_END_A=$( grep "^$LOW_START" $DAT | cut -d " " -f2 )
					L_START_B=$( grep "^$LOW_START" $DAT | cut -d " " -f4 )
					L_END_B=$( grep "^$LOW_START" $DAT | cut -d " " -f5 )
					H_START_A=$( grep " $HIGH_END " $DAT | cut -d " " -f1 )
					H_END_A=$( grep " $HIGH_END " $DAT | cut -d " " -f2 )		
					H_START_B=$( grep " $HIGH_END " $DAT | cut -d " " -f4 )
					echo -e "$LOW_START\t""$L_END_A\t""repeat_region\n""\t\t\trpt_type\tITR\n""$L_START_B\t""$L_END_B\t""repeat_region\n""\t\t\trpt_type\tITR\n""$H_START_A\t""$H_END_A\t""repeat_region\n""\t\t\trpt_type\tITR\n""$H_START_B\t""$HIGH_END\t""repeat_region\n""\t\t\trpt_type\tITR\n" >> ../${DAT%.fasta.2.3.5.80.10.40.500000.10000.dat}.ITR.tbl; 
				fi ; 
			fi ; 
		fi ; 
	done

	# Looking for non-circular contigs without ITRs that have at least 1 virus-specific domain
	echo "$(tput setaf 4) Looking for non-circular contigs without ITRs that have at least 1 virus-specific or plasmid-specific domain $(tput sgr 0)"

	mkdir ../no_end_contigs_with_viral_domain
	for NO_END in $CONTIGS_NON_CIRCULAR ; do 
		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 150 -sequence $NO_END -outseq ${NO_END%.fasta}.AA.fasta
		bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' ${NO_END%.fasta}.AA.fasta > ${NO_END%.fasta}.AA.sorted.fasta ;
		# Taking arguments for "virus specific" database and conducting hmmscan
		### Make an argument to NOT look at non-ITR and non-circular
		if  [[ $virus_domain_db = "-standard" ]] ; then
			hmmscan --tblout ${NO_END%.fasta}.AA.hmmscan.out --cpu 56 -E 1e-8 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/baits_edited_1/virus_specific_baits_edited2 ${NO_END%.fasta}.AA.sorted.fasta
		elif [[ $virus_domain_db = "-with_rdrp_retro" ]]; then
			hmmscan --tblout ${NO_END%.fasta}.AA.hmmscan.out --cpu 56 -E 1e-8 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/baits_edited1_plus_rdrp_retro/baits_edited1_plus_rdrp_retro2_db ${NO_END%.fasta}.AA.sorted.fasta
		elif [[ $virus_domain_db = "-all_common" ]]; then
			hmmscan --tblout ${NO_END%.fasta}.AA.hmmscan.out --cpu 56 -E 1e-8 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/useful_hmms_list3/viral_domains_total_edited2 ${NO_END%.fasta}.AA.sorted.fasta
		elif [[ $virus_domain_db = "-gambia" ]]; then
			hmmscan --tblout ${NO_END%.fasta}.AA.hmmscan.out --cpu 56 -E 1e-8 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/viral_baits_w_gambia/virus_specific_baits_and_gambia1 ${NO_END%.fasta}.AA.sorted.fasta
		else
			echo "$(tput setaf 5) Incorrect argument given for virus_domain_db variable. Try -standard, -with_rdrp_retro, -all_common as arguments. For this run, no contigs with viral domains but without circularity or ITRs will be detected $(tput sgr 0)"
			rm ./*{0..9}.fasta
			break
		fi

		grep -v "^#" ${NO_END%.fasta}.AA.hmmscan.out | sed 's/ \+/	/g' | sort -u -k3,3 > ${NO_END%.fasta}.AA.hmmscan.sort.out
		VIRAL_HALLMARK_COUNT=$( cat ${NO_END%.fasta}.AA.hmmscan.sort.out | wc -l )
		if [ $VIRAL_HALLMARK_COUNT -gt $MINIMUM_DOMAINS ] ; then
			
			echo $NO_END "contains at least $MINIMUM_DOMAINS viral or plasmid domain(s)"

			cut -f3 ${NO_END%.fasta}.AA.hmmscan.sort.out > ${NO_END%.fasta}.rotate.AA.called_hmmscan.txt ; 
			grep -v -f ${NO_END%.fasta}.rotate.AA.called_hmmscan.txt ${NO_END%.fasta}.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.no_hmmscan1.fasta
			echo ">Feature "${NO_END%.fasta}" Table1" > ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.SCAN.tbl

			mv $NO_END ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.fna
			mv ${NO_END%.fasta}.AA.hmmscan.out ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.AA.hmmscan.out
			mv ${NO_END%.fasta}.AA.hmmscan.sort.out ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.AA.hmmscan.sort.out
			#mv ${NO_END%.fasta}.AA.fasta ../no_end_contigs_with_viral_domain/${NO_END%.fasta}.AA.fasta
		else 
			cat $NO_END >> noncircular_non_viral_domains_contigs.fna
			rm $NO_END
		fi
	done
fi


cd ..
echo "$(tput setaf 4) transferring ITR-containing contigs to the main directory for full annotation$(tput sgr 0)"
ITR_ELEMENTS=$( ls ITR_containing_contigs/*.fasta | sed 's/ITR_containing_contigs\///g' )
if [ ! -z "$ITR_ELEMENTS" ] ; then
	for ITR_CONTIG in $ITR_ELEMENTS ; do
		### make .tbl extension of inverted repeat features
		echo $ITR_CONTIG "has ITRs"
		cp ITR_containing_contigs/$ITR_CONTIG $ITR_CONTIG
	done
fi
# Checking whether any circular contigs are >90% identical to any sequence in NCBI nt database using BLASTN.
CIRCLES_AND_ITRS=$( ls | grep ".fasta" )
if [ ! -z "$CIRCLES_AND_ITRS" ] ; then 
	echo "$(tput setaf 5)Checking Novelty of:$(tput sgr 0)"
	for line in $CIRCLES_AND_ITRS ; do
		echo "$(tput setaf 5)"$line"$(tput sgr 0)"
	done
	echo " "
	mkdir circles_of_known_viruses
	mkdir circles_of_chromosomal_elements
	for circle in $CIRCLES_AND_ITRS ; do
		if [[ $handle_knowns = "-quick_knowns" ]] ; then

			echo "starting BLASTN"
			blastn -db /fdb/blastdb/nt -query $circle -evalue 1e-50 -num_threads 56 -outfmt "6 qseqid sseqid stitle pident length qlen" -qcov_hsp_perc 50 -num_alignments 3 -out ${circle%.fasta}.blastn.out ;
			cat ${circle%.fasta}.blastn.out
			if [ -s "${circle%.fasta}.blastn.out" ]; then
				echo ${circle%.fasta}.blastn.out" found"
				sed 's/ /-/g' ${circle%.fasta}.blastn.out | awk '{if ($4 > 90) print}' | awk '{if (($5 / $6) > 0.5) print}' > ${circle%.fasta}.blastn.notnew.out ;
			else
				echo ${circle%.fasta}.blastn.out" not found"
			fi
			if [ -s "${circle%.fasta}.blastn.notnew.out" ] ; then

				echo "$(tput setaf 4)"$circle" is not a novel species (>90% identical to sequence in GenBank nt database).$(tput sgr 0)"
				ktClassifyBLAST -o ${circle%.fasta}.tax_guide.blastn.tab ${circle%.fasta}.blastn.notnew.out
				taxid=$( tail -n1 ${circle%.fasta}.tax_guide.blastn.tab | cut -f2 )
				efetch -db taxonomy -id $taxid -format xml | /data/tiszamj/mike_tisza/xtract.Linux -pattern Taxon -element Lineage > ${circle%.fasta}.tax_guide.blastn.out
				sleep 2s
				if [ !  -z "${circle%.fasta}.tax_guide.blastn.out" ] ; then
					awk '{ print "; "$3 }' ${circle%.fasta}.blastn.notnew.out | sed 's/-/ /g; s/, complete genome//g' >> ${circle%.fasta}.tax_guide.blastn.out
				fi

				if grep -i -q "virus\|viridae\|virales\|Circular-genetic-element\|Circular genetic element\|plasmid" ${circle%.fasta}.tax_guide.blastn.out ; then
					echo $circle "$(tput setaf 4) is closely related to a virus that has already been deposited in GenBank nt and will be annotated with BLASTP only. $(tput sgr 0)"
					cat ${circle%.fasta}.tax_guide.blastn.out
					mv $circle circles_of_known_viruses/${circle%.fasta}.known_species.fna ;
					mv ${circle%.fasta}.tax_guide.blastn.out circles_of_known_viruses/${circle%.fasta}.tax_guide.blastn.out
					mv ${circle%.fasta}.tax_guide.blastn.tab circles_of_known_viruses/${circle%.fasta}.tax_guide.blastn.tab
					mv ${circle%.fasta}.blastn.notnew.out circles_of_known_viruses/${circle%.fasta}.blastn.notnew.out
					mv ${circle%.fasta}.blastn.out circles_of_known_viruses/${circle%.fasta}.blastn.out
				else 
					echo $circle "$(tput setaf 4) is closely related to a chromosomal sequence that has already been deposited in GenBank nt and will be checked for viral and plasmid domains. $(tput sgr 0)"
					cat ${circle%.fasta}.tax_guide.blastn.out
					mv $circle circles_of_chromosomal_elements/${circle%.fasta}.chromosomal.fna
					mv ${circle%.fasta}.tax_guide.blastn.out circles_of_chromosomal_elements/${circle%.fasta}.tax_guide.blastn.out
					mv ${circle%.fasta}.tax_guide.blastn.tab circles_of_chromosomal_elements/${circle%.fasta}.tax_guide.blastn.tab
					mv ${circle%.fasta}.blastn.notnew.out circles_of_chromosomal_elements/${circle%.fasta}.blastn.notnew.out
					mv ${circle%.fasta}.blastn.out circles_of_chromosomal_elements/${circle%.fasta}.blastn.out	
				fi
			else
				echo "$(tput setaf 5)"$circle" appears to be a novel sequence (no close (>90% nucleotide) matches to sequences in nt database).$(tput sgr 0)"
			fi
		elif [[ $handle_knowns = "-do_not_check_knowns" ]] ; then
			echo "$(tput setaf 5) Not checking circular seqs against genbank 'nt' database for close matches $(tput sgr 0)"
		fi
	done


	NOVEL_CIRCLES_AND_ITRS=$( ls | grep "[0-9].fasta" | grep -v "trnascan-se2.txt" )
	echo "$(tput setaf 5)NOVEL_CIRCLES_AND_ITRS variable: $(tput sgr 0)"
	echo $NOVEL_CIRCLES_AND_ITRS
	echo " " 
	# echoing novel circular sequences
	echo " "
fi
if [ -z "$NOVEL_CIRCLES_AND_ITRS" ] ; then
	echo "$(tput setaf 4)No novel circular sequences detected. $(tput sgr 0)" 
	
else
	echo "$(tput setaf 5)Novel circular sequence(s) detected$(tput sgr 0)"
	echo " "
fi

echo "$(tput setaf 5)Automatically annotating novel sequences:$(tput sgr 0)"
for lineg in $NOVEL_CIRCLES_AND_ITRS ; do
	echo "$(tput setaf 5)"$lineg"$(tput sgr 0)"
done
echo " "

# Rotating each sequence to put a non-intragenic start codon as the first basepair of the contig
for nucl_fa in $NOVEL_CIRCLES_AND_ITRS ; do
	if [ -s ITR_containing_contigs/$nucl_fa ] ; then
		cp $nucl_fa ${nucl_fa%.fasta}.rotate.fasta
	else
		echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"

		/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -minsize 240 -find 3 -sequence $nucl_fa -outseq ${nucl_fa%.fasta}.nucl_orfs.fa ; 

		grep ">" ${nucl_fa%.fasta}.nucl_orfs.fa > ${nucl_fa%.fasta}.nucl_orfs.txt
		cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
			start_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
			end_base=$( echo $liner | sed 's/.*- \(.*\)\].*/\1/' )
			length=$(( $start_base-$end_base ))
			abso_length=$( echo $length | sed 's/-//g' )
			if [ $abso_length -gt 239 ]; then
				if [[ "$end_base" -gt "$start_base" ]]; then
					for ((counter_f=(( $start_base + 1 ));counter_f<=(( $end_base + 3 ));counter_f++)); do
						echo "$counter_f" >> ${nucl_fa%.fasta}.bad_starts.txt
						
					done
				elif [[ "$start_base" -gt "$end_base" ]]; then
					for ((counter_r=(( $end_base - 3 ));counter_r<=(( $start_base - 1 ));counter_r++)) ; do
						echo "$counter_r" >> ${nucl_fa%.fasta}.bad_starts.txt
					done
				fi
			fi
		done

		cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
			starter_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
			if grep -q "$starter_base" ${nucl_fa%.fasta}.bad_starts.txt ; then
				continue
			else
				echo $liner >> ${nucl_fa%.fasta}.good_start_orfs.txt
			fi
		done
		if [ -s "${nucl_fa%.fasta}.good_start_orfs.txt" ]; then	
			cut -d " " -f1 ${nucl_fa%.fasta}.good_start_orfs.txt | head -n1 | sed 's/>//g' > ${nucl_fa%.fasta}.starting_orf.txt
			bioawk -c fastx '{ print $name, $seq, length($seq) }' ${nucl_fa%.fasta}.nucl_orfs.fa | grep -f ${nucl_fa%.fasta}.starting_orf.txt | sed '/--/d' | head -n1 | awk '{print ">"$1, $3; print $2}' > ${nucl_fa%.fasta}.starting_orf.1.fa ;
			circlator fixstart --genes_fa ${nucl_fa%.fasta}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fasta}.rotate ;
		else
			cp $nucl_fa ${nucl_fa%.fasta}.no_100AA_ORFs.fasta
		fi
	fi
done 



# Performing BLASTX of each contig against database of viral and plasmid proteins to guess taxonomy
for nucl_fa in $NOVEL_CIRCLES_AND_ITRS ; do
if [ -s "${nucl_fa%.fasta}.rotate.fasta" ]; then
	echo "$(tput setaf 5)Guessing taxonomy for sequence "${nucl_fa%.fasta}.rotate.fasta" by BLASTX against virus and plasmid protein database.$(tput sgr 0)"
	blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -num_threads 56 -num_alignments 1 -db /data/tiszamj/mike_tisza/auto_annotation_pipeline/blast_DBs/virus_refseq_adinto_polinto_clean_plasmid_prot_190911 -query ${nucl_fa%.fasta}.rotate.fasta -out ${nucl_fa%.fasta}.tax_guide.blastx.out ;
	if [ ! -s "${nucl_fa%.fasta}.tax_guide.blastx.out" ]; then
		echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out ;
	elif grep -q "adinto" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
		echo "Adintovirus" >> ${nucl_fa%.fasta}.tax_guide.blastx.out
	elif grep -i -q "polinton-like" ${vd_fa%.fna}.tax_guide.blastx.out ; then
		echo "Polinton-like virus" >> ${vd_fa%.fna}.tax_guide.blastx.out
	else
		echo "$(tput setaf 5)"$nucl_fa" likely represents a novel virus or plasmid. Getting hierarchical taxonomy info.$(tput sgr 0)"
		ktClassifyBLAST -o ${nucl_fa%.fasta}.tax_guide.blastx.tab ${nucl_fa%.fasta}.tax_guide.blastx.out
		taxid=$( tail -n1 ${nucl_fa%.fasta}.tax_guide.blastx.tab | cut -f2 )
		efetch -db taxonomy -id $taxid -format xml | /data/tiszamj/mike_tisza/xtract.Linux -pattern Taxon -element Lineage >> ${nucl_fa%.fasta}.tax_guide.blastx.out
	fi
else
	echo "$(tput setaf 4)"$nucl_fa" could not be rotated. Likely there were no ORFs of at least 100AA.$(tput sgr 0)" 

fi
done


# Extracting ORFs >240bp, (>90bp for inoviruses/plasmids)
for nucl_fa in $NOVEL_CIRCLES_AND_ITRS ; do
if [ -s "${nucl_fa%.fasta}.rotate.fasta" ]; then
	echo "$(tput setaf 5)"$nucl_fa" taxonomy guessed. Continuing to ORF translation...$(tput sgr 0)"

	if [ -s ITR_containing_contigs/$nucl_fa ]; then
		if grep -i -q "Caudovir\|Ackermannvir\|Herellevir\|Corticovir\|Levivir\|Tectivir\|crAss-like virus\|CrAssphage\|Cyanophage\|Microvir\microphage\|Siphoviridae\|Myoviridae\|phage\|Podovir\|Halovir\|sphaerolipovir\|pleolipovir\|plasmid\|Inovir\|Ampullavir\|Bicaudavir\|Fusellovir\|Guttavir\|Ligamenvir\|Plasmavir\|Salterprovir\|Cystovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then

					/data/tiszamj/mike_tisza/PHANOTATE/phanotate.py -f fasta -o ${nucl_fa%.fasta}.phan.fasta ${nucl_fa%.fasta}.rotate.fasta ; 
					/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/transeq -frame 1 -table 11 -sequence ${nucl_fa%.fasta}.phan.fasta -outseq ${nucl_fa%.fasta}.trans.fasta ; COUNTER=0 ; 
					bioawk -c fastx '{print}' ${nucl_fa%.fasta}.trans.fasta | while read LINE ; do 
						START_BASE=$( echo $LINE | sed 's/.*START=\(.*\)\] \[.*/\1/' ) ; 
						ORF_NAME=$( echo $LINE | cut -d " " -f1 | sed 's/\(.*\)\.[0-9].*_1/\1/' ) ; 
						END_BASE=$( echo $LINE | cut -d " " -f1 | sed 's/.*\(\.[0-9].*_1\)/\1/' | sed 's/_1//g; s/\.//g' ) ; 
						ORIG_CONTIG=$( grep ">" $nucl_fa | cut -d " " -f2 ) ; 
						AA_SEQ=$( echo "$LINE" | cut -f2 | sed 's/\*//g' ) ; 
						let COUNTER=COUNTER+1 ; 
						echo ">"${ORF_NAME}"_"${COUNTER} "["$START_BASE" - "$END_BASE"]" $ORIG_CONTIG  ; echo $AA_SEQ ; 
					done > ${nucl_fa%.fasta}.rotate.AA.fasta

			else
				/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 150 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
			fi
	else
		if grep -i -q "Caudovir\|Ackermannvir\|Herellevir\|Corticovir\|Levivir\|Tectivir\|crAss-like virus\|CrAssphage\|Cyanophage\|Microvir\microphage\|Siphoviridae\|Myoviridae\|phage\|Podovir\|Halovir\|sphaerolipovir\|pleolipovir\|plasmid\|Inovir\|Ampullavir\|Bicaudavir\|Fusellovir\|Guttavir\|Ligamenvir\|Plasmavir\|Salterprovir\|Cystovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then

			/data/tiszamj/mike_tisza/PHANOTATE/phanotate.py -f fasta -o ${nucl_fa%.fasta}.phan.fasta ${nucl_fa%.fasta}.rotate.fasta ; 
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/transeq -frame 1 -table 11 -sequence ${nucl_fa%.fasta}.phan.fasta -outseq ${nucl_fa%.fasta}.trans.fasta ; COUNTER=0 ; 
			bioawk -c fastx '{print}' ${nucl_fa%.fasta}.trans.fasta | while read LINE ; do 
				START_BASE=$( echo $LINE | sed 's/.*START=\(.*\)\] \[.*/\1/' ) ; 
				ORF_NAME=$( echo $LINE | cut -d " " -f1 | sed 's/\(.*\)\.[0-9].*_1/\1/' ) ; 
				END_BASE=$( echo $LINE | cut -d " " -f1 | sed 's/.*\(\.[0-9].*_1\)/\1/' | sed 's/_1//g; s/\.//g' ) ; 
				ORIG_CONTIG=$( grep ">" $nucl_fa | cut -d " " -f2 ) ; 
				AA_SEQ=$( echo "$LINE" | cut -f2 | sed 's/\*//g' ) ; 
				let COUNTER=COUNTER+1 ; 
				echo ">"${ORF_NAME}"_"${COUNTER} "["$START_BASE" - "$END_BASE"]" $ORIG_CONTIG  ; echo $AA_SEQ ; 
			done > ${nucl_fa%.fasta}.rotate.AA.fasta
		else
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -find 1 -minsize 150 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
		fi
	fi


	bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' ${nucl_fa%.fasta}.rotate.AA.fasta > ${nucl_fa%.fasta}.rotate.AA.sorted.fasta ;
fi
done

# Conducting HMMER search on curated database

for nucl_fa in $NOVEL_CIRCLES_AND_ITRS ; do
if [ -s "${nucl_fa%.fasta}.rotate.AA.sorted.fasta" ]; then
	echo "$(tput setaf 5)"$nucl_fa" Continuing to HMMSCAN (HMMER) on custom viral conserved protein model database for each ORF...$(tput sgr 0)" 
	hmmscan --tblout ${nucl_fa%.fasta}.rotate.AA.hmmscan.out --cpu 56 -E 1e-6 /data/tiszamj/mike_tisza/auto_annotation_pipeline/refseq_viral_family_proteins/useful_hmms_list3/viral_domains_total_edited2 ${nucl_fa%.fasta}.rotate.AA.sorted.fasta
	echo "$(tput setaf 5)HMMSCAN search of "${nucl_fa%.fasta}.rotate.AA.sorted.fasta" complete.$(tput sgr 0)"
	echo " "
	grep -v "^#" ${nucl_fa%.fasta}.rotate.AA.hmmscan.out | sed 's/ \+/	/g' | sort -u -k3,3 > ${nucl_fa%.fasta}.rotate.AA.hmmscan.sort.out
	cut -f3 ${nucl_fa%.fasta}.rotate.AA.hmmscan.sort.out | awk '{ print $0" " }' > ${nucl_fa%.fasta}.rotate.AA.called_hmmscan.txt ; 
	if [ -s ${nucl_fa%.fasta}.rotate.AA.called_hmmscan.txt ] ; then

		grep -v -f ${nucl_fa%.fasta}.rotate.AA.called_hmmscan.txt ${nucl_fa%.fasta}.rotate.AA.sorted.fasta | grep -A1 ">" | sed '/--/d' > ${nucl_fa%.fasta}.rotate.no_hmmscan.fasta
		echo ">Feature "${nucl_fa%.fasta}" Table1" > ${nucl_fa%.fasta}.SCAN.tbl

		cat ${nucl_fa%.fasta}.rotate.AA.called_hmmscan.txt | while read LINE ; do 
			PROTEIN_INFO=$( grep "$LINE \[" ${nucl_fa%.fasta}.rotate.AA.sorted.fasta ) ; 
			echo $PROTEIN_INFO 
			START_BASEH=$( echo $PROTEIN_INFO | sed 's/.*\[\(.*\) -.*/\1/' ) ; 
			END_BASEH=$( echo $PROTEIN_INFO | sed 's/.*- \(.*\)\].*/\1/' ) ; 
			HMM_INFO=$( grep "$LINE" ${nucl_fa%.fasta}.rotate.AA.hmmscan.out | head -n1 | cut -d " " -f1 | sed 's/-/ /g; s/.*[0-9]\+\///g') ; 
			INFERENCEH=$( echo $HMM_INFO | cut -d " " -f1 ) ; 
			PROTEIN_NAME=$( echo $HMM_INFO | cut -d " " -f2- ) ; 
			echo -e "$START_BASEH\t""$END_BASEH\t""CDS\n""\t\t\tprotein_id\t""lcl|""$LINE\n""\t\t\tproduct\t""$PROTEIN_NAME\n""\t\t\tinference\t""similar to AA sequence:$INFERENCEH" >> ${nucl_fa%.fasta}.SCAN.tbl ; 
		done
	else
		cp ${nucl_fa%.fasta}.rotate.AA.sorted.fasta ${nucl_fa%.fasta}.rotate.no_hmmscan.fasta
	fi
else
	echo "$(tput setaf 4)ORF file for "$nucl_fa" is empty. This circle might have no ORFS over 100AA.$(tput sgr 0)"
fi
done

### chromosomal script goes here
if [ ! -d "circles_of_chromosomal_elements" ]; then
	mkdir circles_of_chromosomal_elements
fi 
CHROMOSOMAL_ELEMENT_LIST=$( ls circles_of_chromosomal_elements/*.fna )
if [ ! -z "$CHROMOSOMAL_ELEMENT_LIST" ] ; then
	echo "$(tput setaf 5)"$nucl_fa" Searching circular sequences that have close identity to chromosomal sequences for viral or plasmid domains which would indicate that they are proviruses or endogenous plasmids $(tput sgr 0)" 
	. /data/tiszamj/mike_tisza/auto_annotation_pipeline/quick_annotation_of_chromosomal_seqs1.sh
else
	echo "$(tput setaf 5) No circles resembling chromosomal elements detected $(tput sgr 0)"
fi


NEW_FASTAS=$( ls | grep "[0-9].fasta" | grep -v ".no_hmmscan.\|.AA.\|.rotate.\|all_non_circular\|.dat" )

echo $NEW_FASTAS
# Conducting RPS-BLAST against CDD on translated ORFs
for nucl_fa in $NEW_FASTAS ; do
if [ -s "${nucl_fa%.fasta}.rotate.no_hmmscan.fasta" ]; then

	echo "$(tput setaf 5)"$nucl_fa" Continuing to RPS-BLAST NCBI CDD domains database for each ORF...$(tput sgr 0)" 
	rpsblast -evalue 1e-4 -num_descriptions 5 -num_threads 56 -line_length 100 -num_alignments 1 -db /data/tiszamj/mike_tisza/auto_annotation_pipeline/cdd_rps_db/Cdd -query ${nucl_fa%.fasta}.rotate.no_hmmscan.fasta -out ${nucl_fa%.fasta}.rotate.AA.rpsblast.out ;
	echo "$(tput setaf 5)RPS-BLAST of "${nucl_fa%.fasta}.rotate.no_hmmscan.fasta" complete.$(tput sgr 0)"
	echo " "
else
	echo "$(tput setaf 4) no ORFs for CD-HIT from "$nucl_fa". all ORFs may have been called with HMMSCAN.$(tput sgr 0)"
	echo " "
fi
done
#rm ${nucl_fa%.fasta}.nucl_orfs.fa ${nucl_fa%.fasta}.rotate.detailed.log ${nucl_fa%.fasta}.rotate.log ${nucl_fa%.fasta}.rotate.promer.promer ${nucl_fa%.fasta}.rotate.promer.contigs_with_ends.fa ${nucl_fa%.fasta}.rotate.prodigal.for_prodigal.fa ${nucl_fa%.fasta}.rotate.prodigal.prodigal.gff;


# Generating tbl file from RPS-BLAST results
echo "$(tput setaf 5) Starting perl script to make tbl from RPS-BLAST output $(tput sgr 0)"

perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/rpsblastreport2tbl_mt_annotation_pipe_biowulf.pl ;

for nucl_fa in $NEW_FASTAS ; do 
if [ -s "${nucl_fa%.fasta}.NT.tbl" ]; then
	echo "$(tput setaf 5)"$nucl_fa" tbl made from RPS-BLAST hits...$(tput sgr 0)"
else
	echo "$(tput setaf 4) RPS-BLAST tbl for "$nucl_fa" not detected.$(tput sgr 0)"
fi
done

# Conducting BLASTP on ORFs unrecognized by RPS-BLAST (nr virus or nr all database)
for feat_tbl1 in *.NT.tbl ; do
grep -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' -B2 $feat_tbl1 | grep "^[0-9]" | awk '{print $1 " - " $2}' > ${feat_tbl1%.NT.tbl}.for_blastp.txt ;
grep -f ${feat_tbl1%.NT.tbl}.for_blastp.txt -A1 ${feat_tbl1%.NT.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta ;
if grep -q "(plasmid)" ${feat_tbl1%.NT.tbl}.tax_guide.blastx.out ; then
	if [ -s "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" ]; then
		echo "$(tput setaf 5)"$nucl_fa" is likely a plasmid... Continuing to BLASTP NCBI nr database for each ORF that had no hits in CDD...$(tput sgr 0)" 
		blastp -evalue 1e-4 -num_descriptions 5 -num_threads 56 -num_alignments 1 -db /fdb/blastdb/nr -query ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta -out ${feat_tbl1%.NT.tbl}.rotate.blastp.out ;
		echo "$(tput setaf 5)BLASTP of "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" complete.$(tput sgr 0)"
	fi
else
	if [ -s "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" ]; then
		echo "$(tput setaf 5)"$nucl_fa" is likely a virus... Continuing to BLASTP NCBI nr database for each ORF that had no hits in CDD...$(tput sgr 0)" 
		blastp -evalue 1e-4 -num_descriptions 5 -num_threads 56 -num_alignments 1 -db /fdb/blastdb/viral -query ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta -out ${feat_tbl1%.NT.tbl}.rotate.blastp.out ;
		echo "$(tput setaf 5)BLASTP of "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" complete.$(tput sgr 0)"
	fi
fi
done

# Generating tbl file from BLASTP results
echo "$(tput setaf 5) Starting perl script to make tbl from BLASTP output $(tput sgr 0)"

perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/blastpreport2tbl_mt_annotation_pipe_biowulf2.pl ;
for feat_tbl1 in *.NT.tbl ; do
if [ -s "${feat_tbl1%.NT.tbl}.BLASTP.tbl" ]; then
	echo "$(tput setaf 5)"${feat_tbl1%.NT.tbl}": tbl made from BLASTP hits. Splitting fasta files for HHsearch...$(tput sgr 0)"
else
	echo "$(tput setaf 4) BLASTP tbl for "${feat_tbl1%.NT.tbl}" not detected.$(tput sgr 0)"
fi
done

# Detecting any tRNAs and making a tbl addenum file
echo "$(tput setaf 5) Looking for tRNAs in dsDNA phage and large dsDNA viruses $(tput sgr 0)"

for GENOME_NAME in $NEW_FASTAS ; do
	tRNAscan-SE -Q -G -o $GENOME_NAME.trnascan-se2.txt ${GENOME_NAME%.fasta}.rotate.fasta
	
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

#	if [ -s $GENOME_NAME.trna.tbl ] ; then
#		cat OTHER_TBL_FILE $GENOME_NAME.trna.tbl > FINAL_TBL_FILE
#	fi

### Put a intermediate tbl combiner here
echo "$(tput setaf 5) combining .tbl files that have been generated so far $(tput sgr 0)"

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
	elif [ -s ${nucl_fa%.fasta}.NT.tbl ] ; then
		cat ${nucl_fa%.fasta}.NT.tbl > ${nucl_fa%.fasta}.int.tbl
		if [ -s ${nucl_fa%.fasta}.trna.tbl ] ; then
			echo -e "\n" >> ${nucl_fa%.fasta}.int.tbl
			cat ${nucl_fa%.fasta}.trna.tbl >> ${nucl_fa%.fasta}.int.tbl
		fi
	elif [ -s ${nucl_fa%.fasta}.SCAN.tbl ] ; then
		cat ${nucl_fa%.fasta}.SCAN.tbl > ${nucl_fa%.fasta}.int.tbl
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


# remove ORFs within ORFs that are 'hypothetical'
echo "$(tput setaf 5) Removing ORFS within ORFs that are 'hypothetical' $(tput sgr 0)"


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
			echo $feat_tbl3
			echo "$loc1_start"
			if [[ "$loc_end" -gt "$loc_start" ]]; then
				gen_len=$(( $loc_end - $loc_start ))

				if [[ "$gen_len" -gt 1000 ]]; then
					continue
				else
					f_end=$(( $loc_end + 1 ))
					f1_end=$( echo " " "$f_end" " ")
					echo "$f1_end"
					if grep -q "$f1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then
						echo "removing hypo"
						echo "$loc1_start" "start, " "$loc_end" "end, " 
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

				echo "$r1_end"
					if grep -q "$r1_end" ${feat_tbl3%.int.tbl}.used_positions.txt ; then
							echo "removing hypo"
							echo "$loc1_start" "start, " "$loc_end" "end, " 
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



# Grabbing ORFs wihout BLASTP hits and separating them into individual files for HHsearch
echo "$(tput setaf 5) Grabbing ORFs wihout BLASTP hits and separating them into individual files for HHsearch $(tput sgr 0)"

for blastp_tbl1 in *.int2.tbl ; do
grep -i -e 'hypothetical protein' -e 'unnamed protein product' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Uncharacterized conserved protein' -e 'unknown' -e 'Uncharacterised protein' -e 'product	gp' -e 'putative phage protein' -B2 $blastp_tbl1 | grep "^[0-9]" | awk '{print $1 " - " $2}' > ${blastp_tbl1%.int2.tbl}.for_hhpred.txt ;
grep -f ${blastp_tbl1%.int2.tbl}.for_hhpred.txt -A1 ${blastp_tbl1%.int2.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${blastp_tbl1%.int2.tbl}.rotate.blast_hypo.fasta ;
csplit -z ${blastp_tbl1%.int2.tbl}.rotate.blast_hypo.fasta '/>/' '{*}' --prefix=${blastp_tbl1%.int2.tbl}.rotate --suffix-format=%02d.for_hhpred.fasta; 
done

# Running HHsearch on remaining ORFs
dark_orf_list=$( ls *.for_hhpred.fasta )

for dark_orf in $dark_orf_list ; do
	if  [[ $HHSUITE_TOOL = "-hhsearch" ]] ; then

		echo "$(tput setaf 5)Running HHsearch on "$dark_orf" now.$(tput sgr 0)"
		hhsearch -i $dark_orf -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pdb70/pdb70 -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pfam_31_db/pfam -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/cdd_db/NCBI_CD/NCBI_CD -o ${dark_orf%.for_hhpred.fasta}.out.hhr -cpu 56 -maxmem 56 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1  ;
		cat ${dark_orf%.for_hhpred.fasta}.out.hhr >> ${dark_orf%.rotate*.for_hhpred.fasta}.rotate.out_all.hhr ;
		rm ${dark_orf%.for_hhpred.fasta}.out.hhr 
		cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
		rm $dark_orf
	elif [[ $HHSUITE_TOOL = "-hhblits" ]] ; then
		echo "$(tput setaf 5)Running HHblits on "$dark_orf" now.$(tput sgr 0)"
		hhblits -i $dark_orf -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pdb70/pdb70 -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pfam_31_db/pfam -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/cdd_db/NCBI_CD/NCBI_CD -o ${dark_orf%.for_hhpred.fasta}.out.hhr -cpu 56 -maxmem 56 -p 80 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1  ;
		cat ${dark_orf%.for_hhpred.fasta}.out.hhr >> ${dark_orf%.rotate*.for_hhpred.fasta}.rotate.out_all.hhr ;
		rm ${dark_orf%.for_hhpred.fasta}.out.hhr 
		cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
		rm $dark_orf
	else
		echo "$(tput setaf 5) Valid option for HHsuite tool (i.e. -hhsearch or -hhblits) was not provided. Skipping step for "$dark_orf" $(tput sgr 0)"
		cat $dark_orf >> ${dark_orf%.rotate*.for_hhpred.fasta}.all_hhpred_queries.AA.fasta
		rm $dark_orf
	fi
done

rm *.rotate.AA.fasta

# Generating tbl file from HHsearch results
echo "$(tput setaf 5) Starting perl script to make tbl from HHsearch output $(tput sgr 0)"

perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/hhpredreport2tbl_mt_annotation_pipe_biowulf1_gjs_edits.pl ;

for HH_tbl1 in *.HH.tbl ; do 
sed 's/OS=.*//g; s/ ;//g; s/UniProtKB:>\([0-9][A-Z].*\)/PDB:\1/g; s/UniProtKB:>tr|.*|\(.\)/UniProtKB:\1/g; s/UniProtKB:>\([a-z].*\)/Scop:\1/g; s/UniProtKB:>\(PF.*\)/PFAM:\1/g; s/ is .*//g; s/ are .*//g' $HH_tbl1 > ${HH_tbl1%.HH.tbl}.HH2.tbl
done

# Combining tbl files from all search results AND fix overlapping ORF module
echo "$(tput setaf 5) Combining tbl files from all search results AND fix overlapping ORF module $(tput sgr 0)"

for feat_tbl4 in *.int2.tbl ; do 
	if [ -s "${feat_tbl4%.int2.tbl}.HH2.tbl" ] && [ -s "$feat_tbl4" ] ; then
		head -n1 $feat_tbl4 > ${feat_tbl4%.int2.tbl}.comb3.tbl
		grep -v -e 'hypothetical protein' -e 'unnamed protein product' -e 'unknown' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -e 'product	gp' $feat_tbl4 | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d; s/TPA_asm: //g' >> ${feat_tbl4%.int2.tbl}.comb3.tbl
		sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g' ${feat_tbl4%.int2.tbl}.HH2.tbl | grep -v ">Feature" | sed '/--/d' >> ${feat_tbl4%.int2.tbl}.comb3.tbl ;
		if [ -s "${feat_tbl4%.int2.tbl}.ITR.tbl" ] ; then
			cat ${feat_tbl4%.int2.tbl}.ITR.tbl >> ${feat_tbl4%.int2.tbl}.comb3.tbl
		fi
	else
		cat $feat_tbl4 | sed '/--/d; s/TPA_asm: //g' > ${feat_tbl4%.int2.tbl}.comb3.tbl
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

rm *tmp.tbl

# .gtf file maker
echo "$(tput setaf 5) Making .gff files for each annotated sequence $(tput sgr 0)"

for feat_tbl2 in *.comb3.tbl ; do
	if [ -s ${feat_tbl2%.comb3.tbl}.gtf ] ; then
		rm ${feat_tbl2%.comb3.tbl}.gtf
	fi
	grep "^[0-9]" -A3 $feat_tbl2 | sed '/--/d' | sed 's/ /_/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n		//g' | while read LINE ; do
		if echo $LINE | grep -q "CDS" ; then
			GENOME=${feat_tbl2%.comb3.tbl}
			FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
			FEAT_START=$( echo $LINE | cut -d " " -f1 )
			FEAT_END=$( echo $LINE | cut -d " " -f2 )
			FEAT_NAME=$( echo $LINE | cut -d " " -f7 )
			FEAT_ATT=$( echo $LINE | cut -d " " -f9 )
			FEAT_ID=$( echo $LINE | cut -d " " -f5 )
		elif echo $LINE | grep -q "repeat_region" ; then
			GENOME=${feat_tbl2%.comb3.tbl}
			FEAT_TYPE=$( echo $LINE | cut -d " " -f3 )
			FEAT_START=$( echo $LINE | cut -d " " -f1 )
			FEAT_END=$( echo $LINE | cut -d " " -f2 )
			FEAT_NAME="ITR"
			FEAT_ATT="ITR"
			FEAT_ID="ITR"	
		fi


		echo -e "$GENOME\t""Cenote-Taker\t""$FEAT_TYPE\t""$FEAT_START\t""$FEAT_END\t"".\t"".\t"".\t""gene_id \"$FEAT_ID\"; gene_name \"$FEAT_NAME\"; gene_inference \"$FEAT_ATT\"" >> ${feat_tbl2%.comb3.tbl}.gtf
	done
done

# Making directory for sequin generation
if [ ! -d "sequin_directory" ]; then
	mkdir sequin_directory
fi

# Getting info for virus nomenclature and divergence 
echo "$(tput setaf 5) Getting info for virus nomenclature and divergence $(tput sgr 0)"

for feat_tbl2 in *.comb3.tbl ; do 
	file_core=${feat_tbl2%.comb3.tbl}
	echo $file_core
	file_numbers=$( echo ${file_core: -3} | sed 's/[a-z]//g' | sed 's/[A-Z]//g' )
	echo $file_numbers
	tax_info=${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out
	echo $tax_info
	if grep -q "Anellovir" $tax_info ; then
		vir_name=Anelloviridae ;
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
	elif grep -q "Cyanophage" $tax_info ; then
		vir_name=Cyanophage ;
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

	elif grep -q "No homologues found" $tax_info ; then
		if  [ -s ITR_containing_contigs/${feat_tbl2%.comb3.tbl}.fasta ] ; then
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
	elif grep -q "Satellite" $tax_info ; then
		vir_name="Satellite virus" ;
	elif grep -q "unclassified ssDNA bacterial virus" $tax_info ; then
		vir_name="unclassified ssDNA bacterial virus" ;
	elif grep -q "phage" $tax_info ; then
		vir_name="Phage" ;
	elif grep -q "plasmid" $tax_info ; then
		vir_name="metagenomic plasmid" ;
	elif grep -q "Bacteria" $tax_info ; then
		vir_name="Phage" ;
	elif grep -q "virus" $tax_info ; then
		vir_name="Virus" ;
	else
		if  [ -s ITR_containing_contigs/${feat_tbl2%.comb3.tbl}.fasta ] ; then
			vir_name="unclassified element" ;
		else
			vir_name="Circular genetic element" ;
		fi
	fi
	echo $vir_name ;
#	feature_head=$( echo ">Feature" $3"-associated" $vir_name $file_numbers" Table1" )
	fsa_head=$( echo $vir_name " sp." )
	tax_guess=$( tail -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out ) ; 
	perc_id=$( head -n1 ${feat_tbl2%.comb3.tbl}.tax_guide.blastx.out | sed 's/ /-/g' | awk '{FS="\t"; OFS="\t"} {print $2" "$3}' | sed 's/-/ /g' ) ;
	rand_id=$( echo $RANDOM | tr '[0-9]' '[a-zA-Z]' | cut -c 1-2 )

	# Editing and transferring tbl file and fasta (fsa) files to sequin directory
	echo "$(tput setaf 5) Editing and transferring tbl file and fasta (fsa) files to sequin directory $(tput sgr 0)"

	if [ -s ${feat_tbl2%.comb3.tbl}.phan.fasta ]; then
		echo "$(tput setaf 5)tbl file made from: "$feat_tbl2" will be used for sqn generation$(tput sgr 0)" ; 
		cp $feat_tbl2 sequin_directory/${feat_tbl2%.comb3.tbl}.tbl ; 
		if [ -s ITR_containing_contigs/${feat_tbl2%.comb3.tbl}.fasta ]; then
			if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then
				bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var " ; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"

			else
				bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
			fi

		else
			if [ -s circles_of_chromosomal_elements/${feat_tbl2%.comb3.tbl}.provirus.fna ] ; then
				CELL_CHROM=$( cat ${feat_tbl2%.comb3.tbl}.blastn.notnew.out | head -n1 | cut -f3 )
				if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v chrom_info="$CELL_CHROM" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var " ; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [note=highly similar to sequence "chrom_info "; please manually check if this is a transposon especially if there is an annotated reverse transcriptase] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
					echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"

				else
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v chrom_info="$CELL_CHROM" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [note=highly similar to sequence "chrom_info "; please manually check if this is a transposon especially if there is an annotated reverse transcriptase] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				fi
			else
				if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then				
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
					echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"
				else
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				fi
			fi
		fi
	else
		echo "$(tput setaf 5)tbl file made from: "$feat_tbl2" will be used for sqn generation$(tput sgr 0)" ; 
		cp $feat_tbl2 sequin_directory/${feat_tbl2%.comb3.tbl}.tbl ; 
		if [ -s ITR_containing_contigs/${feat_tbl2%.comb3.tbl}.fasta ]; then
			if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then
				bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var " ; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"

			else
				bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
			fi

		else
			if [ -s circles_of_chromosomal_elements/${feat_tbl2%.comb3.tbl}.provirus.fna ] ; then
				CELL_CHROM=$( cat ${feat_tbl2%.comb3.tbl}.blastn.notnew.out | head -n1 | cut -f3 )
				if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v chrom_info="$CELL_CHROM" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var " ; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [note=highly similar to sequence "chrom_info "; please manually check if this is a transposon especially if there is an annotated reverse transcriptase] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=11]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
					echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"

				else
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v chrom_info="$CELL_CHROM" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [note=highly similar to sequence "chrom_info "; please manually check if this is a transposon especially if there is an annotated reverse transcriptase] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				fi
			else
				if [ -z "${feat_tbl2%.comb3.tbl}.rotate.AA.hmmscan.sort.out" ] ; then				
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "; WARNING: no viral/plasmid-specific domains were detected. This is probably not a true mobile genetic element.] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
					echo "$(tput setaf 3) WARNING: no viral/plasmid-specific domains were detected in "${feat_tbl2%.comb3.tbl}". This is probably not a true mobile genetic element. Scrutinize carefully. $(tput sgr 0)"
				else
					bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -v molecule_var="$MOLECULE_TYPE" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername " ct" rand_var number_var "] [moltype=genomic "molecule_var"][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "] [gcode=1]" ; print $seq }' ${feat_tbl2%.comb3.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.comb3.tbl}.fsa ;	
				fi
			fi
		fi
	fi
done


#making cmt file for assembly data
for nucl_fa in $NEW_FASTAS ; do
	echo "making cmt file!!!!!!"
	input_contig_name=$( head -n1 $nucl_fa | cut -d " " -f 1 | sed 's/|.*//g; s/>//g' ) 
	echo $input_contig_name
	COVERAGE=$( grep "$input_contig_name	" reads_to_all_contigs_over${length_cutoff}nt.coverage.txt | cut -f2 )
	echo $COVERAGE
	echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Assembly Method	" $ASSEMBLER >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Genome Coverage	"$COVERAGE"x" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Sequencing Technology	Illumina" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Annotation Pipeline	Cenote-Taker2" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "URL	https://github.com/mtisza1/Cenote-Taker2" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
done


for fsa_file in sequin_directory/*.fsa ; do
	echo "editing fsa file!!!!!!"
	fsa_name2=$( echo ${fsa_file#sequin_directory/} ) ; 
	fsa_name3=$( echo ${fsa_name2%.fsa} | sed 's/.PLASMID//g' )
	seq_name1=$( head -n1 $fsa_name3.fasta | sed 's/>//g; s/|.*//g' | cut -d " " -f2 )
	sed " 1 s/note= closest relative/note= $seq_name1 ; closest relative/" $fsa_file > $fsa_file.temp
	mv $fsa_file.temp $fsa_file
done


# Running sequin to generate sqn, gbf, and val files for each genome
if [[ $DATA_SOURCE = "-tpa_assembly" ]] ;then
	/data/tiszamj/python/linux64.tbl2asn -V vb -j "[keyword=TPA:assembly]" -t $base_directory/$template_file -X C -p sequin_directory/ ;
else
	/data/tiszamj/python/linux64.tbl2asn -V vb -t $base_directory/$template_file -X C -p sequin_directory/ ;
fi

# Script for annotating complete circular viruses of known species
KNOWN_VIRAL_CIRCLES=$( ls circles_of_known_viruses/ | grep ".fna" )
if [ ! -z "$KNOWN_VIRAL_CIRCLES" ] ;then
	echo "$(tput setaf 3) Starting annotation of known circular viruses $(tput sgr 0)"
	. /data/tiszamj/mike_tisza/auto_annotation_pipeline/quick_annotation_of_known_seqs1.sh
fi

# script for annotating no_end contigs with viral domains
LIST_OF_VIRAL_DOMAIN_CONTIGS=$( ls no_end_contigs_with_viral_domain/  | grep ".fna" )
if [ ! -z "$LIST_OF_VIRAL_DOMAIN_CONTIGS" ] ;then
	echo "$(tput setaf 3) Starting annotation of known contigs with viral domains but are neither circular nor have ITRs $(tput sgr 0)"

	. /data/tiszamj/mike_tisza/auto_annotation_pipeline/annotation_of_viral_domain_contigs.sh
fi

# script for sketching contigs without detectable virus domains or circularity or ITRs with RPS-BLAST
if  [[ $handle_nonviral = "-sketch_all" ]] ; then
	echo "$(tput setaf 3) Sketching contigs without detectable virus domains or circularity or ITRs with RPS-BLAST $(tput sgr 0)"
	. /data/tiszamj/mike_tisza/auto_annotation_pipeline/sketching_nonviral_contigs.sh
elif [[ $handle_nonviral = "-no_sketch_domainless" ]]; then
	echo "$(tput setaf 3) The -no_sketch option was used, contigs without detectable virus domains or circularity or ITRs will not be sketched $(tput sgr 0)"
elif [[ $handle_nonviral = "-full_annotate_domainless" ]]; then
	echo "$(tput setaf 3) Fully annotating all contigs without detectable virus domains, circularity, or ITRs with RPS-BLAST $(tput sgr 0)"
	. /data/tiszamj/mike_tisza/auto_annotation_pipeline/full_annotate_of_putative_nonviral_contigs1.sh
else
	echo "$(tput setaf 3) No recognized sketch option was used, contigs without detectable virus domains or circularity or ITRs will not be sketched $(tput sgr 0)"
fi

cd $base_directory/$run_title

# making a summary table (.tsv) for all putatively viral contigs/genomes

echo "$(tput setaf 3) Making a summary table of all viral contigs, if any. $(tput sgr 0)"

echo -e "Isolation source""\t""Completeness""\t""Cenote-taker contig name""\t""Length""\t""Element Name""\t""Topology""\t""Common Viral Domains""\t""BLAST hit for Taxonomy" > ${run_title}.tsv
for i in sequin_directory/*.fsa ; do
	if [ ! -z "$i" ] ;then
		site=$( head -n1 $i | sed -e 's/.*isolation_source=\(.*\)\] \[isolate.*/\1/' )
		echo $site
		df_num=$( head -n1 $i | sed -e 's/>\(.*[0-9].*\)\ \[note= .*/\1/' )
		echo $df_num
		topologyq=$( head -n1 $i | sed -e 's/.*topology=\(.*\)\] \[Bioproject.*/\1/' )
		echo $topologyq
		tax_call=$( head -n1 $i | sed -e 's/.*organism=\(.*\)\] \[moltype=.*/\1/' )
		echo $tax_call
		blast_call=$( head -n1 $i | sed -e 's/.*closest relative: \(.*\)\] \[organism=.*/\1/' )
		echo $blast_call
		length=$( bioawk -c fastx '{ print length}' $i )
		echo length
		j=${i%.fsa} ; 
		DOMAIN_COUNT=$( cat ${j#sequin_directory/}.rotate.AA.hmmscan.sort.out | wc -l )
		#title=$( cat $site | awk '{print $1}' )
		echo -e $site "\t""Complete genome (putative)""\t" $df_num "\t" $length "\t" $tax_call "\t" $topologyq "\t" $DOMAIN_COUNT "\t" $blast_call >> ${run_title}.tsv
	fi
done

for i in circles_of_known_viruses/sequin_directory/*.fsa ; do
	if [ ! -z "$i" ] ;then
		site=$( head -n1 $i | sed -e 's/.*isolation_source=\(.*\)\] \[isolate.*/\1/' )
		echo $site
		df_num=$( head -n1 $i | sed -e 's/>\(.*[0-9].*\)\ \[note= .*/\1/' )
		echo $df_num
		topologyq=$( head -n1 $i | sed -e 's/.*topology=\(.*\)\] \[Bioproject.*/\1/' )
		echo $topologyq
		tax_call=$( head -n1 $i | sed -e 's/.*organism=\(.*\)\] \[moltype=.*/\1/' )
		echo $tax_call
		blast_call1=$( head -n1 $i | sed -e 's/.*closest relative: \(.*\)\] \[organism=.*/\1/' )
		blast_call2=$( echo "CLOSE SIMILARITY TO:" $blast_call1 )
		length=$( bioawk -c fastx '{ print length}' $i )
		echo $length
		DOMAIN_COUNT="Not Determined"
		echo -e $site "\t""Complete genome (putative)""\t" $df_num "\t" $length "\t" $tax_call "\t" $topologyq "\t" $DOMAIN_COUNT "\t" $blast_call2 >> ${run_title}.tsv
	fi
done

for i in no_end_contigs_with_viral_domain/sequin_directory/*.fsa ; do
	if [ ! -z "$i" ] ;then
		site=$( head -n1 $i | sed -e 's/.*isolation_source=\(.*\)\] \[isolate.*/\1/' )
		echo $site
		df_num=$( head -n1 $i | sed -e 's/>\(.*[0-9].*\)\ \[note= this.*/\1/' )
		echo $df_num
		topologyq=$( head -n1 $i | sed -e 's/.*topology=\(.*\)\] \[Bioproject.*/\1/' )
		echo $topologyq
		### INSERT OPTION FOR BLASTN RESULT
		tax_call=$( head -n1 $i | sed -e 's/.*organism=\(.*\)\] \[moltype=.*/\1/' )
		echo $tax_call
		if grep -q "highly similar to sequence" $i ; then
			blast_call1=$( head -n1 $i | sed -e 's/note=highly similar to sequence \(.*\); please manually check if this is a transposon especially if there is an annotated reverse transcriptase.*/\1/' )
			blast_call2=$( echo "CLOSE SIMILARITY TO:" $blast_call1" check if this is a RT transposon" )
		else
			blast_call2=$( head -n1 $i | sed -e 's/.*closest relative: \(.*\)\] \[organism=.*/\1/' )
		fi
		echo $blast_call
		length=$( bioawk -c fastx '{ print length}' $i )
		echo $length
		j=${i%.fsa} ; 
		DOMAIN_COUNT=$( cat no_end_contigs_with_viral_domain/${j#no_end_contigs_with_viral_domain/sequin_directory/}.AA.hmmscan.sort.out no_end_contigs_with_viral_domain/${j#no_end_contigs_with_viral_domain/sequin_directory/}.AA.hmmscan2.sort.out | sort -k1,1 | wc -l )
		echo -e $site "\t""Partial genome (putative)""\t" $df_num "\t" $length "\t" $tax_call "\t" $topologyq "\t" $DOMAIN_COUNT "\t" $blast_call2 >> ${run_title}.tsv
	fi
done

echo "$(tput setaf 3) Summary file made: ${run_title}.tsv $(tput sgr 0)"


rm *.all_start_stop.txt *.bad_starts.txt *.comb.tbl *.comb2.tbl *.good_start_orfs.txt *.hypo_start_stop.txt *.nucl_orfs.fa *.remove_hypo.txt *.log *.promer.contigs_with_ends.fa *.promer.promer *.out.hhr *.starting_orf.txt *.out.hhr
rm noncircular_contigs/*.AA.fasta noncircular_contigs/*.AA.sorted.fasta noncircular_contigs/*.out noncircular_contigs/*.dat

echo " "
echo "$(tput setaf 3) "$run_title" $(tput sgr 0)"
echo "$(tput setaf 3) >>>>>>CENOTE TAKER HAS FINISHED TAKING CENOTES<<<<<< $(tput sgr 0)"



