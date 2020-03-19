# Cenote-Taker2
Cenote-Taker2 is a pipeline for divergent virus discovery and annotation. See schematic.
The code is currently functional.
An ulterior motive for creating and distributing Cenote-Taker2 is to facilitate annotation and deposition of viral genomes into GenBank where they can be used by the scientific public. Therefore, I hope you consider depositing the submittable outputs (.sqn) after reviewing them. I am not affiliated with GenBank. See the [Cenote-Taker2 wiki](https://github.com/mtisza1/Cenote-Taker2/wiki) for useful information on using the pipeline (e.g. expected outputs) and screeds on myriad topics.

![alt text](https://github.com/mtisza1/Cenote-Taker2/blob/master/cenote-taker_logo.png)

# Install Using Conda

```
**The Five Commandments
**1. Know where your default conda environment install space is.
**2. Ensure you have space for all 130GB of files.
**3. Don't install without checking conda version first.
**4. Only install on an HPC running on Linux, unless your personal computer has sick specs.
**5. Don't let your computer fall asleep during the install
```

Likely, this will only work in Linux. 
Using a HPC with at least 16 CPUs and 16g of dedicated memory is strongly recommended. (Annotation of a few selected genomes can be done with less memory/CPU) 
I usually use 32 CPUs and 32 GB of memory for medium and large metagenomes. More resources would be helpful for extra-large metagenomes.
```diff
- ALERT *** Because Cenote-Taker2 needs large high-quality 
- sequence databases to work correctly, installation will take ~2 hours 
- AND require about 130GB of storage space. 
- Also part of the install requires 4 CPUs to work. 
- Therefore, you may need to be in an interactive job on an HPC
```
1. Change to the directory you'd like to be the parent to the install directory
2. Ensure Conda is installed and working. Use version 4.8.2 or better.
```
conda -V
```
3. Download the script in the \'install_scripts\' directory of this github repo into your current directory. (i.e. cenote_install1.sh). (remove any older versions of cenote_install1.sh first, if applicable)
```
wget  https://raw.githubusercontent.com/mtisza1/Cenote-Taker2/master/install_scripts/cenote_install1.sh
```
4. Run the install script. Give exactly one argument in the script: 'default' - OR - a path to the desired conda environment setup directory. 
The conda environment itself requires 32GB of space mostly due to the krona taxonomy database. Some HPC users have installed their Conda in their /home directory (or equivalent) which typically has little space. The other 100GB consists of sequence databases and will be put within the current directory
```diff
- ALERT *** Because Cenote-Taker2 needs large high-quality 
- sequence databases to work correctly, running this script will take ~2 hours 
- AND require about 130GB of storage space. 
- Also part of the install requires 4 CPUs to work. 
- Therefore, you may need to be in an interactive job on an HPC
```

```
If there is enough space in your default conda environment directory:
bash cenote_install1.sh default

Otherwise specify an absolute path to a directory with >32GB of storage:
bash cenote_install1.sh /path/to/better/directory
```


![alt text](https://github.com/mtisza1/Cenote-Taker2/blob/master/cenote-taker_figure_200318.png)

# Running Cenote-Taker2
Cenote-Taker2 currently runs in a python wrapper. 
1. Activate the Conda environment.
```
Default:
conda activate cenote-taker2_env

Or if you've put your conda environment in a custom location:
conda activate /path/to/better/directory/cenote-taker2_env
```
2. Run the python script (see options below).
```
python /path/to/Cenote-Taker2/run_cenote-taker2.0.1.py
```


```
usage: run_cenote-taker2.0.1.py [-h] 
                          --contigs ORIGINAL_CONTIGS 
                          --run_title RUN_TITLE 
                          --template_file TEMPLATE_FILE
                          --prune_prophage PROPHAGE 
                          --mem MEM 
                          --cpu CPU
                                [--reads1 F_READS] 
                                [--reads2 R_READS]
                                [--minimum_length_circular CIRC_LENGTH_CUTOFF]
                                [--minimum_length_linear LINEAR_LENGTH_CUTOFF]
                                [--virus_domain_db VIRUS_DOMAIN_DB]
                                [--lin_minimum_hallmark_genes LIN_MINIMUM_DOMAINS]
                                [--circ_minimum_hallmark_genes CIRC_MINIMUM_DOMAINS]
                                [--known_strains HANDLE_KNOWNS]
                                [--blastn_db BLASTN_DB]
                                [--enforce_start_codon ENFORCE_START_CODON]
                                [--handle_contigs_without_hallmark HANDLE_NONVIRAL]
                                [--hhsuite_tool HHSUITE_TOOL]
                                [--isolation_source ISOLATION_SOURCE]
                                [--Environmental_sample ENVIRONMENTAL_SAMPLE]
                                [--collection_date COLLECTION_DATE]
                                [--metagenome_type METAGENOME_TYPE]
                                [--srr_number SRR_NUMBER]
                                [--srx_number SRX_NUMBER]
                                [--biosample BIOSAMPLE]
                                [--bioproject BIOPROJECT]
                                [--assembler ASSEMBLER]
                                [--molecule_type MOLECULE_TYPE]
                                [--data_source DATA_SOURCE]
                                [--filter_out_plasmids FILTER_PLASMIDS]
                                [--scratch_directory SCRATCH_DIR]
                                [--blastp BLASTP]

Cenote-Taker2 is a pipeline for virus discovery and thorough annotation of
viral contigs and genomes.

optional arguments:
  -h, --help            show this help message and exit

 REQUIRED ARGUMENTS for Cenote-Taker2 :
  --contigs ORIGINAL_CONTIGS
                        Contig file with .fasta extension in fasta format - OR
                        - assembly graph with .fastg extension. Each header
                        must be unique before the first space character
  --run_title RUN_TITLE
                        Name of this run. A directory of this name will be
                        created. Must be unique from older runs or older run
                        will be renamed. Must be less than 18 characters,
                        using ONLY letters, numbers and underscores (_)
  --template_file TEMPLATE_FILE
                        Template file with some metadata. Takes a couple
                        minutes to generate: https://submit.ncbi.nlm.nih.gov/g
                        enbank/template/submission/
  --prune_prophage PROPHAGE
                        True or False. Attempt to identify and remove flanking
                        chromosomal regions from non-circular contigs with
                        viral hallmarks (True is highly recommended for
                        sequenced material not enriched for viruses. Virus
                        enriched samples probably should be False (you might
                        check with ViromeQC). Also, please use False if
                        --lin_minimum_hallmark_genes is set to 0)
  --mem MEM             example: 56 Gigabytes of memory available for Cenote-
                        Taker2. Typically, 16 to 32 should be used. Lower
                        memory will work in theory, but could extend the
                        length of the run
  --cpu CPU             Example: 32 Number of CPUs available for Cenote-
                        Taker2. Typically, 32 CPUs should be used. For large
                        datasets, increased performance can be seen up to 120
                        CPUs. Fewer than 16 CPUs will work in theory, but
                        could extend the length of the run

 OPTIONAL ARGUMENTS for Cenote-Taker2. Most of which are important to consider!!! GenBank typically only accepts genome submission with ample metadata. See https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#ModifiersPage for more information on GenBank metadata fields:
  --reads1 F_READS      Default: no_reads ILLUMINA READS ONLY: First Read file
                        in paired read set - OR - read file in unpaired read
                        set - OR - read file of interleaved reads. Used for
                        coverage depth determination.
  --reads2 R_READS      Default: no_reads ILLUMINA READS ONLY: Second Read
                        file in paired read set. Disregard if not using paired
                        reads. Used for coverage depth determination.
  --minimum_length_circular CIRC_LENGTH_CUTOFF
                        Default: 1000 Minimum length of contigs to be checked
                        for circularity. Bare minimun is 1000 nts
  --minimum_length_linear LINEAR_LENGTH_CUTOFF
                        Default: 1000 Minimum length of non-circualr contigs
                        to be checked for viral hallmark genes.
  --virus_domain_db VIRUS_DOMAIN_DB
                        default: standard 'standard' database: mostly DNA
                        virus hallmark genes (i.e. genes with known function
                        and exclusively found in viral genomes) with very low
                        false discovery rate. 'rna_virus' database: For RNA
                        sequencing data only. Includes RdRp and capsid genes
                        of RNA viruses. Moderate false discovery rate due to
                        structural similarity between RdRp genes and e.g.
                        transposon-encoded RT genes
  --lin_minimum_hallmark_genes LIN_MINIMUM_DOMAINS
                        Default: 1 Number of detected viral hallmark genes on
                        a non-circular contig to be considered viral and
                        recieve full annotation. WARNING: Only choose '0' if
                        you have prefiltered the contig file to only contain
                        putative viral contigs (using another method such as
                        VirSorter or DeepVirFinder), or you are very confident
                        you have physically enriched for virus particles very
                        well (you might check with ViromeQC). Otherwise, the
                        duration of the run will be extended many many times
                        over, largely annotating non-viral contigs, which is
                        not what Cenote-Taker2 is meant for. For unenriched
                        samples, '2' might be more suitable, yielding a false
                        positive rate near 0.
  --circ_minimum_hallmark_genes CIRC_MINIMUM_DOMAINS
                        Default:1 Number of detected viral hallmark genes on a
                        circular contig to be considered viral and recieve
                        full annotation. For samples physically enriched for
                        virus particles, '0' can be used, but please treat
                        circular contigs without known viral domains
                        cautiously. For unenriched samples, '1' might be more
                        suitable.
  --known_strains HANDLE_KNOWNS
                        Default: do_not_check_knowns -> do not check if
                        putatively viral contigs are highly related to known
                        sequences (via MEGABLAST). 'blast_knowns': REQUIRES '
                        --blastn_db' option to function correctly.
  --blastn_db BLASTN_DB
                        Default: none Set a database if using '--
                        known_strains' option. Specify BLAST-formatted
                        nucleotide datase. Probably, use only GenBank 'nt'
                        database downloaded from ftp://ftp.ncbi.nlm.nih.gov/
                        or another GenBank formatted .fasta file to make
                        databse
  --enforce_start_codon ENFORCE_START_CODON
                        Default: True For final genome maps, require ORFs to
                        be initiated by a typical start codon? GenBank
                        submissions containing ORFs without start codons can
                        be rejected. However, if True, important but
                        incomplete genes could be culled from the final
                        output. This is relevant mainly to contigs of
                        incomplete genomes
  --handle_contigs_without_hallmark HANDLE_NONVIRAL
                        Default: no_sketch_domainless What do you want to do
                        with contigs that do not have detectable viral
                        hallmark features? 'no_sketch_domainless': do nothing,
                        report sequences in file. 'sketch_all': annotate ORFs
                        with RPSBLAST/CDD and tRNA scan only (Could still add
                        substantial time to run, especially without phyical or
                        computational viral enrichment).
  --hhsuite_tool HHSUITE_TOOL
                        default: hhblits hhblits will query PDB, pfam, and CDD
                        to annotate ORFs escaping identification via upstream
                        methods. 'hhsearch': hhsearch, a more sensitive tool,
                        will query PDB, pfam, and CDD to annotate ORFs
                        escaping identification via upstream methods.
                        (WARNING: hhsearch takes much, much longer than
                        hhblits and can extend the duration of the run many
                        times over. Do not use on large input contig files).
                        'no_hhsuite_tool': forgoes annotation of ORFs with
                        hhsuite. Fastest way to complete a run.
  --isolation_source ISOLATION_SOURCE
                        Default: unknown Describes the local geographical
                        source of the organism from which the sequence was
                        derived
  --Environmental_sample ENVIRONMENTAL_SAMPLE
                        Default: False True or False, Identifies sequence
                        derived by direct molecular isolation from an
                        unidentified organism
  --collection_date COLLECTION_DATE
                        Default: unknown Date of collection. this format:
                        01-Jan-2019, i.e. DD-Mmm-YYYY
  --metagenome_type METAGENOME_TYPE
                        Default: unknown a.k.a. metagenome_source
  --srr_number SRR_NUMBER
                        Default: unknown For read data on SRA, run number,
                        usually beginning with 'SRR' or 'ERR'
  --srx_number SRX_NUMBER
                        Default: unknown For read data on SRA, experiment
                        number, usually beginning with 'SRX' or 'ERX'
  --biosample BIOSAMPLE
                        Default: unknown For read data on SRA, sample number,
                        usually beginning with 'SAMN' or 'SAMEA' or 'SRS'
  --bioproject BIOPROJECT
                        Default: unknown For read data on SRA, project number,
                        usually beginning with 'PRJNA' or 'PRJEB'
  --assembler ASSEMBLER
                        Default: unknown_assembler Assembler used to generate
                        contigs, if applicable. Specify version of assembler
                        software, if possible.
  --molecule_type MOLECULE_TYPE
                        Default: DNA viable options are DNA - OR - RNA
  --data_source DATA_SOURCE
                        default: original original data is not taken from
                        other researchers' public or private database.
                        'tpa_assembly': data is taken from other researchers'
                        public or private database. Please be sure to specify
                        SRA metadata.
  --filter_out_plasmids FILTER_PLASMIDS
                        Default: True True - OR - False. If True, hallmark
                        genes of plasmids will not count toward the minimum
                        hallmark gene parameters. If False, hallmark genes of
                        plasmids will count. Plasmid hallmark gene set is not
                        necessarily comprehensive at this time.
  --scratch_directory SCRATCH_DIR
                        Default: none When running many instances of Cenote-
                        Taker2, it seems to run more quickly if you copy the
                        hhsuite databases to a scratch space temporarily. Use
                        this argument to set a scratch directory that the
                        databases will be copied to (at least 100GB of scratch
                        space are required for copying the databases)
  --blastp BLASTP       Do not use this argument as of now.
```

![alt text](https://github.com/mtisza1/Cenote-Taker2/blob/master/cenote-taker2_directory_tree2.png)
