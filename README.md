# Cenote-Taker 2
Cenote-Taker 2 is a dual function bioinformatics tool. On the one hand, Cenote-Taker 2 can discover/predict virus sequences from any kind of genome or metagenomic assembly. On the other hand, virus sequences/genomes (perhaps predicted by another tool?) can be annotated with a variety of sequences features, genes, and taxonomy. Either the discovery or the the annotation module can be used independently.
```diff
+ The code is currently functional. Feel free to use Cenote-Taker 2 at will.
+ Major update on May 6th 2022: Version 2.1.5
+ Cenote-Taker 2.1.5 has an easier, more reliable installation and database downloads. Some packages that have given many users issues have been replaced. Taxonomy is more flexible. See release notes.
```

If you just want to discover/predict virus sequences and get a report on those sequences, use [Cenote Unlimited Breadsticks](https://github.com/mtisza1/Cenote_Unlimited_Breadsticks), also provided in the Cenote-Taker 2 repo.

If you just want to annotate your virus sequences and make genome maps, run Cenote-Taker 2 using `-am True`.

An ulterior motive for creating and distributing Cenote-Taker 2 is to facilitate annotation and deposition of viral genomes into GenBank where they can be used by the scientific public.  Therefore, I hope you consider depositing the submittable outputs (.sqn) after reviewing them. I am not affiliated with GenBank. 

See "Use Cases" below, and read the [Cenote-Taker 2 wiki](https://github.com/mtisza1/Cenote-Taker2/wiki) for useful information on using the pipeline (e.g. expected outputs) and screeds on myriad topics.
Using a HPC with at least 16 CPUs and 16g of dedicated memory is recommended for most runs. (Annotation of a few selected genomes or virus discovery on smaller databases can be done with less memory/CPU in a reasonable amount of time). 

To update from `v2.1.3`  (note that biopython and bedtools are now required): 
```
conda activate cenote-taker2_env
conda install -c phanotate=1.5.0 seqkit=2.2.0 hhsuite=3.3.0
cd Cenote-Taker2
git pull
#Then update the BLAST database (see instructions below).
```

Update to HMM databases (hallmark genes) occurred on June 16th, 2021. Update to the BLAST (taxonomy) database occurred on May 6th, 2022. See instructions below to update your database.

Read the manuscript in [Virus Evolution](https://academic.oup.com/ve/article/7/1/veaa100/6055568)

If you cannot or do not want to install and run this on the command line, `Cenote-Taker 2 v 2.1.3` is freely available to run with point-and-clink interface on the [CyVerse Discovery Environment](https://cyverse.org/discovery-environment).


![alt text](../master/cenote-taker_logo.png)

# Install Using Conda

```

** Databases will require between 8GB (most basic) and 75GB (all the optional databases) of storage.
** Don't install without checking conda version first.
** Install on machine running on Linux (with a reasonably new OS), or a Mac (see instructions).
```

If you just want a lightweight (7GB), faster, **NON-ANNOTATING** virus discovery tool, use [Cenote Unlimited Breadsticks](https://github.com/mtisza1/Cenote_Unlimited_Breadsticks). The `Unlimited Breadsticks` module is included in the Cenote-Taker 2 repo, so no need to install it if you already have `Cenote-Taker 2` (you may need to update from older versions `Cenote-Taker2`) 

```diff
- ALERT *** If you choose to install all optional databases for HHsuite, 
- installation will take about 2 hours due to slow download speeds for pdb70
- AND require about 75GB of storage space. 
```
1. Change to the directory you'd like to be the parent to the install directory

2. Ensure Conda is installed and working (required for installation and execution of `Cenote-Taker 2`). Use version 4.10 or better.
Note: instructions for installing Conda are probably specific to your university's/organization's requirements, so it is always best to ask your IT professional or HPC administrator. Generally, you will want to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) in your data directory.
```
conda -V
```
3. Clone the `Cenote-Taker 2` github repo.

```
git clone https://github.com/mtisza1/Cenote-Taker2.git
```
4. Install the conda environment (phanotate and hhsuite don't play nice with the .yml file, so they need special commands)

```
conda env create --file cenote-taker2_env.yml
# follow conda prompts to allow install

conda activate cenote-taker2_env

pip install phanotate

conda install -c conda-forge -c bioconda hhsuite last=1282
```

5. Change to the `Cenote-Taker2` repo directory OR a different location where you want the databases to be stored. (NOTE: if you install the databases in a custom location you will need to specify this directory each time you run the tool) Download the databases.

```
conda activate cenote-taker2_env
cd Cenote-Taker2

# with all the options (75GB)
python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True --hhPDB True

# substantially smaller but with some hhsuite DBs (20GB)
python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True

# only the required DBs, No hhsuite (8GB)
python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True
``` 


# Bioconda installation

* THIS HAS NOT BEEN UPDATED RECENTLY. BIOCONDA VERSION NOT RECOMMENDED AT THE MOMENT *

A user has packaged Cenote-Taker 2 in Bioconda for use by their institute.  However, installation can be done by anyone using their package with a few commands. All the above alerts, requirements, and warnings still apply. This will also require a user to have 32GB of storage in their default conda environment directory.

Commands:
```
conda create -n cenote-taker2 -c hcc -c conda-forge -c bioconda -c defaults cenote-taker2=2020.04.01

conda activate cenote-taker2

download-db.sh
```

The Krona database directory will then need to be manually downloaded and set up. This should work:
```
CT2_DIR=$PWD
KRONA_DIR=$( which python | sed 's/bin\/python/opt\/krona/g' )
cd ${KRONA_DIR}
sh updateTaxonomy.sh
cd ${KRONA_DIR}
sh updateAccessions.sh
cd ${CT2_DIR}
```
Discussion:
[LINK](https://github.com/mtisza1/Cenote-Taker2/issues/1#issuecomment-608510204)

## Updating databases

As of now, the HMM database has been updated from the original (update on June 16th, 2021), and the BLAST database has been updated (May 6th, 2022). This update should only take a minute or two. Here's how you update (modify if your conda environment is different than below example):
update Cenote-Taker 2 (change to main repo directory):
`git pull`

load your conda environment:
`conda activate cenote-taker2_env`

run the update script:
`python update_ct2_databases.py --hmm True --protein True`


## Schematic
![alt text](../master/CT2_schematic_redo1.png)

# Running Cenote-Taker 2
Cenote-Taker 2 currently runs in a python wrapper. 
1. Activate the Conda environment.

Check environments:

`conda info --envs`

Default:

`conda activate cenote-taker2_env`

Or if you've put your conda environment in a custom location:

`conda activate /path/to/better/directory/cenote-taker2_env`

2. Run the python script to get the help menu (see options below).

```
# quick help menu
python /path/to/Cenote-Taker2/run_cenote-taker2.py

# full help menu
python /path/to/Cenote-Taker2/run_cenote-taker2.py -h
```

3. Run some contigs. For example:
```
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_CONTIGS.fasta -r my_contigs1_ct -m 32 -t 32 -p true -db virion

#Or, if you want to save a log of the run, add  "2>&1 | tee output.log" to the end of the command):

python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_CONTIGS.fasta -r my_contigs1_ct -m 32 -t 32 -p true -db virion 2>&1 | tee output.log
```


### Use Case Suggestions/Settings
#### *Annotation*

If you just want to annotate your pre-selected virus sequences and make genome maps, run Cenote-Taker 2 using `-am True`.

Example:
```
# clip and wrap circular sequences
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_VIRUSES.fasta -r viruses_am_ct -m 32 -t 32 -p False -am True

# do not wrap circular sequences, but label DTR regions
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_VIRUSES.fasta -r viruses_am_ct -m 32 -t 32 -p False -am True --wrap True
```

For very divergent genomes, setting `-hh hhsearch` will marginally improve number of genes that are annotated. This setting increasese the run time quite a bit. On the other hand, setting `-hh none` will skip the time consuming hhblits step. With this, you'll still get pretty good genome maps, and might be most appropriate for very large virus genome databases, or for runs where you just want to do a quick check.

#### *Discovery*

**Virus-like particle (VLP) prep assembly:**

`-p False -db standard`

You might apply a size cutoff for linear contigs as well, e.g. ` --minimum_length_linear 3000` OR `--minimum_length_linear 5000`. Changing length minima does not affect false positive rates, but short linear contigs may not be useful, depending on your goals.

Example:
```
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_VLP_ASSEMBLY.fasta -r my_VLP1_ct -m 32 -t 32 -p False -db standard --minimum_length_linear 3000
```

**Whole genome shotgun (WGS) metagenomic assembly:**

`-p True -db virion --minimum_length_linear 3000 --lin_minimum_hallmark_genes 2`

While you should definitely ***definitely*** prune virus sequences from WGS datasets, [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) also does a very good job (I'm still formally comparing these approaches) and you could use `--prune_prophage False` on a metagenome assembly and feed the unpruned contigs from Unlimited Breadsticks into `checkv end_to_end` if you prefer. My suggestion is to prune with `Cenote-Taker 2`, then run `CheckV`.

Example with prune:
```
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_WGS_ASSEMBLY.fasta -r my_WGS1_ct -m 32 -t 32 -p True -db virion --minimum_length_linear 3000 --lin_minimum_hallmark_genes 2
```


**Bacterial reference genome or MAG**

`-p True -db virion --minimum_length_linear 3000 --lin_minimum_hallmark_genes 2`

Using `--lin_minimum_hallmark_genes 1 -db virion` with WGS or bacterial genome data will (in my experience) yield very few sequences that appear to be false positives, however, there are lots of "degraded" prophage sequences in these sequencing sets, i.e. some/most genes of the phage have been lost. That said, sequence with just 1 hallmark gene is neither a guarantee of a degraded phage (especially in the case of ssDNA viruses) nor is 2+ hallmark a guarantee of of a complete phage.

Example:
```
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_BACTERIAL_GENOME.fasta -r my_genome1_ct -m 32 -t 32 -p True -db virion --minimum_length_linear 3000 --lin_minimum_hallmark_genes 2
```


**RNAseq assembly of any kind (if you only want RNA viruses)**

`-p False -db rna_virus`

If you also want DNA virus transcripts, or if your data is mixed RNA/DNA sequencing, you might do a run with `-db rna_virus`, then, from this run, take the file "other_contigs/non_viral_domains_contigs.fna" and use it as input for another run with `-db virion`.

Example:
```
python /path/to/Cenote-Taker2/run_cenote-taker2.py -c MY_METATRANSCRIPTOME.fasta -r my_metatrans1_ct -m 32 -t 32 -p False -db rna_virus
```

All arguments:
```
usage: run_cenote-taker2.py [-h] 
                            -c ORIGINAL_CONTIGS 
                            -r RUN_TITLE 
                            -p PROPHAGE
                            -m MEM 
                            -t CPU 

                            [-am ANNOTATION_MODE]
                            [--template_file TEMPLATE_FILE] 
                            [--reads1 F_READS]
                            [--reads2 R_READS]
                            [--minimum_length_circular CIRC_LENGTH_CUTOFF]
                            [--minimum_length_linear LINEAR_LENGTH_CUTOFF]
                            [-db VIRUS_DOMAIN_DB]
                            [--lin_minimum_hallmark_genes LIN_MINIMUM_DOMAINS]
                            [--circ_minimum_hallmark_genes CIRC_MINIMUM_DOMAINS]
                            [--known_strains HANDLE_KNOWNS]
                            [--blastn_db BLASTN_DB]
                            [--enforce_start_codon ENFORCE_START_CODON]
                            [-hh HHSUITE_TOOL] 
                            [--crispr_file CRISPR_FILE]
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
                            [--orf-within-orf ORF_WITHIN]
                            [--cenote-dbs CENOTE_DBS] [--wrap WRAP]
                            [--hallmark_taxonomy HALLMARK_TAX]


Cenote-Taker 2 is a pipeline for virus discovery and thorough annotation of viral contigs and genomes. 
Visit https://github.com/mtisza1/Cenote-Taker2#use-case-suggestionssettings for suggestions about how to 
run different data types and https://github.com/mtisza1/Cenote-Taker2/wiki to read more. Version 2.1.5

optional arguments:

  -h, --help            show this help message and exit



 REQUIRED ARGUMENTS for Cenote-Taker2 :

  -c ORIGINAL_CONTIGS, --contigs ORIGINAL_CONTIGS

                        Contig file with .fasta extension in fasta format - OR

                        - assembly graph with .fastg extension. Each header

                        must be unique before the first space character

  -r RUN_TITLE, --run_title RUN_TITLE

                        Name of this run. A directory of this name will be

                        created. Must be unique from older runs or older run

                        will be renamed. Must be less than 18 characters,

                        using ONLY letters, numbers and underscores (_)

  -p PROPHAGE, --prune_prophage PROPHAGE

                        True or False. Attempt to identify and remove flanking

                        chromosomal regions from non-circular contigs with

                        viral hallmarks (True is highly recommended for

                        sequenced material not enriched for viruses. Virus

                        enriched samples probably should be False (you might

                        check with ViromeQC). Also, please use False if

                        --lin_minimum_hallmark_genes is set to 0)

  -m MEM, --mem MEM     example: 56 -- Gigabytes of memory available for

                        Cenote-Taker2. Typically, 16 to 32 should be used.

                        Lower memory will work in theory, but could extend the

                        length of the run

  -t CPU, --cpu CPU     Example: 32 -- Number of CPUs available for Cenote-

                        Taker2. Approximately 32 CPUs should be used

                        moderately sized metagenomic assemblies. For large

                        datasets, increased performance can be seen up to 120

                        CPUs. Fewer than 16 CPUs will work in theory, but

                        could extend the length of the run. See GitHub repo

                        for suggestions.



 OPTIONAL ARGUMENTS for Cenote-Taker2. Most of which are important to consider!!! 
 GenBank typically only accepts genome submission with ample metadata. 
 See https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#ModifiersPage for more information on GenBank metadata fields:

  -am ANNOTATION_MODE, --annotation_mode ANNOTATION_MODE

                        Default: False -- Annotate sequences only (skip

                        discovery). Only use if you believe each provided

                        sequence is viral

  --template_file TEMPLATE_FILE

                        Template file with some metadata. Real one required

                        for GenBank submission. Takes a couple minutes to

                        generate: https://submit.ncbi.nlm.nih.gov/genbank/temp

                        late/submission/

  --reads1 F_READS      Default: no_reads -- ILLUMINA READS ONLY: First Read

                        file in paired read set - OR - read file in unpaired

                        read set - OR - read file of interleaved reads. Used

                        for coverage depth determination.

  --reads2 R_READS      Default: no_reads -- ILLUMINA READS ONLY: Second Read

                        file in paired read set. Disregard if not using paired

                        reads. Used for coverage depth determination.

  --minimum_length_circular CIRC_LENGTH_CUTOFF

                        Default: 1000 -- Minimum length of contigs to be

                        checked for circularity. Bare minimun is 1000 nts

  --minimum_length_linear LINEAR_LENGTH_CUTOFF

                        Default: 1000 -- Minimum length of non-circualr

                        contigs to be checked for viral hallmark genes.

  -db VIRUS_DOMAIN_DB, --virus_domain_db VIRUS_DOMAIN_DB

                        default: virion -- 'standard' database: all virus (DNA

                        and RNA) hallmark genes (i.e. genes with known

                        function as virion structural, packaging, replication,

                        or maturation proteins specifically encoded by virus

                        genomes) with low false discovery rate. 'virion'

                        database: subset of 'standard', hallmark genes

                        encoding virion structural proteins, packaging

                        proteins, or capsid maturation proteins (DNA and RNA

                        genomes) with LOWEST false discovery rate. 'rna_virus'

                        database: For RNA virus hallmarks only. Includes RdRp

                        and capsid genes of RNA viruses. Low false discovery

                        rate.

  --lin_minimum_hallmark_genes LIN_MINIMUM_DOMAINS

                        Default: 1 -- Number of detected viral hallmark genes

                        on a non-circular contig to be considered viral and

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

                        Default:1 -- Number of detected viral hallmark genes

                        on a circular contig to be considered viral and

                        recieve full annotation. For samples physically

                        enriched for virus particles, '0' can be used, but

                        please treat circular contigs without known viral

                        domains cautiously. For unenriched samples, '1' might

                        be more suitable.

  --known_strains HANDLE_KNOWNS

                        Default: do_not_check_knowns -- do not check if

                        putatively viral contigs are highly related to known

                        sequences (via MEGABLAST). 'blast_knowns': REQUIRES '

                        --blastn_db' option to function correctly.

  --blastn_db BLASTN_DB

                        Default: none -- Set a database if using '--

                        known_strains' option. Specify BLAST-formatted

                        nucleotide datase. Probably, use only GenBank 'nt'

                        database downloaded from ftp://ftp.ncbi.nlm.nih.gov/

                        or another GenBank formatted .fasta file to make

                        databse

  --enforce_start_codon ENFORCE_START_CODON

                        Default: False -- For final genome maps, require ORFs

                        to be initiated by a typical start codon? GenBank

                        submissions containing ORFs without start codons can

                        be rejected. However, if True, important but

                        incomplete genes could be culled from the final

                        output. This is relevant mainly to contigs of

                        incomplete genomes

  -hh HHSUITE_TOOL, --hhsuite_tool HHSUITE_TOOL

                        default: hhblits -- hhblits will query PDB, pfam, and

                        CDD to annotate ORFs escaping identification via

                        upstream methods. 'hhsearch': hhsearch, a more

                        sensitive tool, will query PDB, pfam, and CDD to

                        annotate ORFs escaping identification via upstream

                        methods. (WARNING: hhsearch takes much, much longer

                        than hhblits and can extend the duration of the run

                        many times over. Do not use on large input contig

                        files). 'no_hhsuite_tool': forgoes annotation of ORFs

                        with hhsuite. Fastest way to complete a run.

  --crispr_file CRISPR_FILE

                        Tab-separated file with CRISPR hits in the following

                        format: CONTIG_NAME HOST_NAME NUMBER_OF_MATCHES. You

                        could use this tool:

                        https://github.com/edzuf/CrisprOpenDB. Then reformat

                        for Cenote-Taker 2

  --isolation_source ISOLATION_SOURCE

                        Default: unknown -- Describes the local geographical

                        source of the organism from which the sequence was

                        derived

  --Environmental_sample ENVIRONMENTAL_SAMPLE

                        Default: False -- True or False, Identifies sequence

                        derived by direct molecular isolation from an

                        unidentified organism

  --collection_date COLLECTION_DATE

                        Default: unknown -- Date of collection. this format:

                        01-Jan-2019, i.e. DD-Mmm-YYYY

  --metagenome_type METAGENOME_TYPE

                        Default: unknown -- a.k.a. metagenome_source

  --srr_number SRR_NUMBER

                        Default: unknown -- For read data on SRA, run number,

                        usually beginning with 'SRR' or 'ERR'

  --srx_number SRX_NUMBER

                        Default: unknown -- For read data on SRA, experiment

                        number, usually beginning with 'SRX' or 'ERX'

  --biosample BIOSAMPLE

                        Default: unknown -- For read data on SRA, sample

                        number, usually beginning with 'SAMN' or 'SAMEA' or

                        'SRS'

  --bioproject BIOPROJECT

                        Default: unknown -- For read data on SRA, project

                        number, usually beginning with 'PRJNA' or 'PRJEB'

  --assembler ASSEMBLER

                        Default: unknown_assembler -- Assembler used to

                        generate contigs, if applicable. Specify version of

                        assembler software, if possible.

  --molecule_type MOLECULE_TYPE

                        Default: DNA -- viable options are DNA - OR - RNA

  --data_source DATA_SOURCE

                        default: original -- original data is not taken from

                        other researchers' public or private database.

                        'tpa_assembly': data is taken from other researchers'

                        public or private database. Please be sure to specify

                        SRA metadata.

  --filter_out_plasmids FILTER_PLASMIDS

                        Default: True -- True - OR - False. If True, hallmark

                        genes of plasmids will not count toward the minimum

                        hallmark gene parameters. If False, hallmark genes of

                        plasmids will count. Plasmid hallmark gene set is not

                        necessarily comprehensive at this time.

  --scratch_directory SCRATCH_DIR

                        Default: none -- When running many instances of

                        Cenote-Taker2, it seems to run more quickly if you

                        copy the hhsuite databases to a scratch space

                        temporarily. Use this argument to set a scratch

                        directory that the databases will be copied to (at

                        least 100GB of scratch space are required for copying

                        the databases)

  --blastp BLASTP       Do not use this argument as of now.

  --orf-within-orf ORF_WITHIN

                        Default: False -- Remove called ORFs without HMMSCAN

                        or RPS-BLAST hits that begin and end within other

                        ORFs? True or False

  --cenote-dbs CENOTE_DBS

                        Default: cenote_script_path -- If you downloaded and

                        setup the databases in a non-standard location,

                        specify path

  --wrap WRAP           Default: True -- Wrap/rotate DTR/circular contigs so

                        the start codon of an ORF is the first nucleotide in

                        the contig/genome

  --hallmark_taxonomy HALLMARK_TAX

                        Default: False -- Get hierarchical taxonomy

                        information for all hallmark genes? This report

                        (*.hallmarks.taxonomy.out) is not considered in the

                        final taxonomy call.
```
## Directory Tree
![Directory Tree Image](../master/cenote-taker2_directory_tree2.png)


## Citation
Michael J Tisza, Anna K Belford, Guillermo Domínguez-Huerta, Benjamin Bolduc, Christopher B Buck, Cenote-Taker 2 democratizes virus discovery and sequence annotation, Virus Evolution, Volume 7, Issue 1, January 2021, veaa100, https://doi.org/10.1093/ve/veaa100

