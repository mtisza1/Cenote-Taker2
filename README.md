# Cenote-Taker2
Cenote-Taker2

![alt text](https://github.com/mtisza1/Cenote-Taker2/blob/master/cenote-taker2_schematic_190920.png)


| Cenote-Taker2 arguments required | That's right, it needs all 21 arguments |
| ------------------------------- | -------- |
| Note	| For all arguments that contain a space (e.g. "gut metagenome") you must use quotes around the phrase |
| 1	| **File of contigs** , MUST end in .fasta or .fastg. NEVER EVER use raw reads, only contigs. |
| 2	| **Forward (1) reads**. Full path to reads. If supplying multiple read sets use format (with quotes): "/data/reads/reads_A_1.fastq, /data/reads/reads_B_1.fastq, /data/reads/reads_C_1.fastq". If reads are not available, you must use "-no_reads" argument here! |
| 3	| **Reverse (2) reads**. Full path to reads. If supplying multiple read sets use format (with quotes): "/data/reads/reads_A_2.fastq, /data/reads/reads_B_2.fastq, /data/reads/reads_C_2.fastq" If reads are not available, you must use "-no_reads" argument here! ALTERNATIVELY: if using unpaired Illumina reads use "-unpaired". If using Nanopore reads list them in the Forward reads spaces and use "-nanopore" here. If using Pacbio reads list them in the Forward reads spaces and use "-pacbio" here. |
| 4	| **Run title**. Must be less than 18 characters, using letters, numbers and underscores (_) |
| 5	| **Isolation source**, e.g.: "cow rumen", see https://www.ncbi.nlm.nih.gov/books/NBK53701/ |
| 6	| **Collection date**, using this format: 01-Jan-2019 (DD-MMM-YYYY) |
| 7	| **Metagenome type**, e.g.: "environmental metagenome", see https://www.insdc.org/documents/feature-table under metagenome_source |
| 8	| **SRR number**, usually begins with "SRR"* |
| 9	| **SRX number**, usually begins with "SRX"* |
| 10	| **Biosample**, usually begins with "SAMN"* |
| 11	| **Bioproject**, usually begins with "PRJNA"* |
| 12	| **Template file**, should be in the working directory. Can be created here https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ |
| 13	| **How to handle putative non-viral contigs.** -sketch_all option will annotate ORFs with RPS-BLAST and find tRNAs. -full_annotate_domainless option will do a full annotation of all contigs regardless of whether a viral domain was detected USE WITH EXTREME CAUTION AND NEVER USE WITH METAGENOMES THAT WERE NOT WELL-ENRICHED FOR VIRAL SEQUENCES. -no_sketch option will not annotate contigs with viral domains (recommended for beginners) |
| 14	| **Length Cutoff**, minimum length of contig for Cenote-Taker to consider, I usually do 1000 for virus enriched datasets. |
| 15	| **HMM database for detecting viral hallmark genes in non-circular non-ITR-containing contigs**. -standard contains only HMMs that I believe won't ping any non-viral contigs. -with_rdrp_retro is -standard plus retrotranscriptase and RNA-dependent RNA polymerase; Good for RNA viromes, but will pull out retro-transposons and the like. |
| 16	| **Minimum domains for non-circular non-ITR-containing contigs**, minimum number of ORFs with hallmark viral domains for further consideration. I recommend 2. However, 1 can be used if you are knowledgable and careful. |
| 17	| **Handle known, complete viral genomes in dataset.** -quick_knowns will only use BLASTP to annotate ORFs, -do_not_check_knowns will not screen viral contigs for close nucleotide identity to viral genomes in GenBank |
| 18	| **Assembler** (include version when known), e.g. "SPAdes v3.1.2" |
| 19	| **Molecule type**, use "DNA" or "RNA" |
| 20	| **HHSuite tool**, -hhsearch will use Hhsearch (TAKES A LONG TIME, BUT IS MOST SENSITIVE), -hhblits is several-fold faster, but less sensitive, -no_hhsuite will skip this process (fastest way to finish, but often forgoes annotation of important ORFs) |
| 21	| **Data source**, use -original for sequences derived from reads that your lab depositied. Use -tpa_assembly for sequences from another lab's reads |
	
	* If you don't currently have SRA information for your samples, just use SRRXXX, SRXXXX, etc. However, you almost certainly will not be able to deposit these outputs in genbank with this info. You can add the information later using Cenote-Fixer.
