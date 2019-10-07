Analyse scripts for A new adaptive index to predict patterns of adaptive genetic variation from SNP-environment associations
============================================================================================================================

Pierre-Edouard Guerin

2017

This folder contains all the scripts to calculate statistics (nucleotide diversity, Tajima's D...) on the genome from SNP data

## Prerequisite
You must install the following softwares :

### VCFTOOLS
[http://vcftools.sourceforge.net/]
  apt-get install vcftools
###TABIX & BGZIP
[https://github.com/samtools/htslib/releases/tag/1.4.1]
  apt-get install tabix
### R Version 3.2.3
[https://cran.r-project.org/]
  install.packages("ggplot2")
  install.packages("PopGenome")
### Python 2.7.12
[https://www.python.org/]

## Data Files
The included data files are :

* `data/Positions.14409.txt`             : List of ID|position|scaffold|chromosome of 14409 SNPs.
* `data/Data.950.sauvages.txt`           : List of names of the 950 indivuals of interest.
* `data/NoPool.14409.csv`                : Table of genotypes of all the indivuals for the 14409 SNPs.
* `data/Noms.marqueurs.LFMM.gINLAnd.csv` : list of SNP IDs and results from Ginland & LFMM methods
* `bon_exemple.vcf`                      : Template of VCF format file used to create VCF files.


## Scripts code sources
scripts used to calculate statistics on the genome from SNP data

### BASH scripts
* `vcf4PopGenome_protocole.sh` : Creates VCF.GZ files using [TABIX], [BGZIP] and [VCFTOOLS] for each chromosome into *mes_vcf/* folder. Uncompressed VCF files will be saved into a new *mes_vcf_save/* folder.
* `fabrique_outlier.sh`        : Generates a list of outliers SNPs positions for each chromosome into *mes_outliers/* folder.
* `add_ID_to_tables.sh` : Creates tables with ID of SNPs into *tables/avec_id* folder from *tables/* results.

### Python scripts
* `get_col.py`          : Select columns of a CSV file according to a list of colunm's names.
* `convert_data2vcf.py` : Creates VCF files for each chromosome of each SNP with genotype of each indivuals.
* `get_id_snp.py`       : get position|chromosome information in a table of SNP and find his ID in a VCF file.

### R scripts
* `fabrique_outlier.R` : From `data/Noms.marqueurs.LFMM.gINLAnd.csv`, it provides a list of outliers SNPs IDs (`nom_84_outliers.txt`)
* `generate_fig_tab.R` : Generates sliding windows genome statistics into *figures/* & *tables/* folders using *mes_vfs/* & *mes_outliers/* data.
* `analysis_tables.R`  : Basic statistical analysis on `tables/avec_id/all_stats_propre.csv`

## Workflows

### Calculate statistics (nucleotide diversity, Tajima's D...) on outliers and non-outliers SNPs and analysis

1. Creates VCF files of all the SNP and selected individuals
  * Run `vcf4PopGenome_protocole.sh` to create VCF.GZ files into *mes_vcf/* folder and VCF files into *mes_vcf_save/*, using *data/* files and `bon_exemple.vcf`

2. Get outliers SNPs and their positions
  * Run `fabrique_outlier.R` to create `nom_84_outliers.txt`, the list of outliers SNPs IDs
  * Create a new *mes_outliers/* folder
  * Run `fabrique_outlier.sh` to create a list of outliers SNPs positions for each chromosome into *mes_outliers/* folder

3. Generates figures and table of genome sliding windows statistics
  * Run `generate_fig_tab.R` to generate plot .PDF figures into *figures/* folder and .CSV tables into *tables/* for each chromosome

4. Analyse statistics on tables
  * Create an empty folder *tables/avec_id/*
  * Run `add_ID_to_tables.sh` to generate a merged table `tables/avec_id/all_stats_propre.csv` with SNPs IDs and statistics of all the chromosomes from *tables/* to *tables/avec_id/*
  * Run `analysis_tables.R` to have outliers SNP and non-outliers SNP statistics, then do some basic statistical tests.


