#!/usr/bin/bash

fix_type=$1
template_file=$2

mkdir cenote_fixer_sequin

if  [[ $fix_type = "-non_rotated" ]] ; then
	cp *.fsa cenote_fixer_sequin/
	cp *.cmt cenote_fixer_sequin/
	cp *.tbl cenote_fixer_sequin/
	$HOME/tbl2asn -V vb -t $template_file -X C -p cenote_fixer_sequin/ ;
elif [[ $fix_type = "-rotated" ]] ; then
	cp *.fsa cenote_fixer_sequin/
	cp *.cmt cenote_fixer_sequin/

	for GENBANK_FILE in *.gb ; do
		counter_f=1
		GENOME_NAME=$( echo $GENBANK_FILE | sed 's/.gb//g; s/_annotations.*//g' )

		echo ">Feature "$GENOME_NAME" Table1" > cenote_fixer_sequin/${GENOME_NAME}.tbl

		grep "     CDS" $GENBANK_FILE | while read LINE ; do
			CDS_FEATURE=$( grep -A3 "$LINE" $GENBANK_FILE )
			PRODUCT=$( echo $CDS_FEATURE | sed 's/.*product="\(.*\)"/\1/' )
			countFound=$( echo $CDS_FEATURE | awk '$0 ~ /inference/ { print }' )
			if [ ! -z "${countFound}" ] ; then
				INFERENCE=$( echo $CDS_FEATURE | sed 's/.*inference="\(.*\)" .*/\1/' )
			else
				INFERENCE=$( echo $CDS_FEATURE | sed 's/.*note="\(.*\)" .*/\1/' )
			fi
			if echo $LINE | grep -q "complement" ; then
				START_BASE=$( echo $CDS_FEATURE | sed 's/.*complement(\(.*\)\.\.\(.*\)).*/\2/' )
				END_BASE=$( echo $CDS_FEATURE | sed 's/.*complement(\(.*\)\.\.\(.*\)).*/\1/' )

			else
				START_BASE=$( echo $CDS_FEATURE | cut -d "." -f1 | cut -d " " -f2 )
				END_BASE=$( echo $CDS_FEATURE | cut -d "." -f3 | cut -d " " -f1 )
			fi
			counter_f=$(( $counter_f + 1 ))
#			echo "feature "$CDS_FEATURE
#			echo "product "$PRODUCT
#			echo "inference "$INFERENCE
#			echo "starting base "$START_BASE
#			echo "ending base "$END_BASE
#			echo "---------"
			if echo $LINE | grep -q "inference" ; then
				echo -e "$START_BASE\t""$END_BASE\t""CDS\n""\t\t\tprotein_id\t""lcl|""$GENOME_NAME""_$counter_f""\n""\t\t\tproduct\t""$PRODUCT\n""\t\t\tinference\t""$INFERENCE" >> cenote_fixer_sequin/${GENOME_NAME}.tbl
			else
				echo -e "$START_BASE\t""$END_BASE\t""CDS\n""\t\t\tprotein_id\t""lcl|""$GENOME_NAME""_$counter_f""\n""\t\t\tproduct\t""$PRODUCT\n""\t\t\tnote\t""$INFERENCE" >> cenote_fixer_sequin/${GENOME_NAME}.tbl
			fi


		done
		$HOME/tbl2asn -V vb -t $template_file -X C -p cenote_fixer_sequin/ ;

	done
fi



