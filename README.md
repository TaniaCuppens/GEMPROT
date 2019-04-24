# GEMPROT
Visualization of the impact on the protein of genetic variants found on each haplotype

## Introduction
GEMPROT (Genetic Mutation to Protein Translation) is a tool to visualize the changes induced on the protein by the different variants found within a gene in an individual. Starting from a phased vcf, it translates the two haplotypes of each individual into the two corresponding protein sequences, allowing to vizualize if variants are in cis or trans. 
GEMPROT proposes two different modes: the first mode outputs both sequences for each individuals and provides a frequency summary of all haplotypes. The second mode is intended for larger populations and shows difference of genes haplotypes repartition.

## Table of content

-	[Dependencies](#Dependencies)
-	[Configuration](#Configuration)
-	[Running GEMPROT](#Running-GEMPROT)
-	[GEMPROT Input](#GEMPROT-Input)
-	[GEMPROT Output](#GEMPROT-Output)
-	[RunExamples](#RunExamples)
-	[Contact Information](#Contact-Information)
-	[Publication](#Publication)

## Dependencies

GEMPROT is written in Perl v5.22.1 and requires the following packages before running:
- perl with "libxml2" package 
`sudo apt-get install libxml-libxml-perl`
-	[Vcftools](https://sourceforge.net/projects/vcftools/files/) 
-	[Samtools](https://sourceforge.net/projects/samtools/files/samtools/0.1.19/)  
-	R with "optparse" package
`install.packages("optparse")`

GEMPROT run with a reference genome (fasta) and its index (.fai), with two CCDS files (Consensus Coding Sequence) and with clinvar data file. These files already exist in the download directory but if you want to use others files please change their path in the configuration file. 
You can download in the FTP sites (be careful about taking the data on the same reference genome): 
GRCh37 :
- ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15/CCDS.current.txt
- ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15/CCDS2Sequence.current.txt
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20180429.vcf.gz  

GRCh38 :
- ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
- ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2Sequence.current.txt
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20180429.vcf.gz


GEMPROT needs an internet access to import domain information via [Pfam](https://pfam.xfam.org/) (Protein families database).

## Configuration


GEMPROT needs some path information.  

conf/conf.pm : 

Reference
############

$genome_reference= "/FILE_PATH/*.fasta";

$CCDS_file="/FILE_PATH/CCDS.current.txt";

$CCDS_NM_file="/FILE_PATH/CCDS2Sequence.current.txt";

$clin_var="/FILE_PATH/clinvar_*.vcf.gz";

Tools
#########

$samtools="/DIR_PATH/samtools";

$vcftools_bin_dir="/DIR_PATH/bin/";

## Running GEMPROT

```perl Script/haplotype_translation_final.pl -h

 2 modes :
        -indiv
    MANDATORY:
                --phased-vcf    : phased vcf
                --output-dir    : specify output directory
                --sample | --sample-file        : sample or sample file (one sample name by line)
                --gene | --gene-file    : official gene symbol or official gene symbol with one ccds in tab delimited file (ex: TP53    CCDS11118.1)
                --assembly (GRCh37 / hg19 or GRCh38 / hg38)
        -pop
    MANDATORY:
                --phased-vcf    : phased vcf
                --output-dir    : specify output directory
                --sample-file   : sample name with his location ; tab delimited file (ex: HG00096       EUR)
                --gene | --gene-file    : official gene symbol or official gene symbol with one ccds in tab delimited file (ex: TP53    CCDS11118.1)
                


    OPTIONAL:
            --loc-file      : location file
            --fasta : if you want protein fasta file for each haplotype and reference
            --synonymous    : if you want to see synonymous SNP
            --domain : if you know protein domain
            -h      : show this message and quit
```
## Input files

The input file is required to be in vcf vcf.gz format. This file can contain one or more individuals and either biallelic or multiallelic sites. You can either provide a reference sequence in the fasta format or a file containing a list of gene name and then choose your coding sequence code by its NM_ reference.

## Output files

Indiv mode
Genename/Visual/
1.	Genename_sample_protein_plot.pdf : visualization of each haplotype in gene with his known domain
2.	Sample_genename_result.html : resume of the analysis

Pop mode
Genename/Visual/
1.	Genename_hap_protein_plot.pdf : View of each haplotype in the population
2.	Result_ analysis.html : Resume of the analysis

With `-fasta` option:

Genename/Fasta/
1.	Sample.sequencegene.ref.prot.fa : reference protein sequence
2.	Sample.sequencegene.hap1.prot.fa : haplotype 1 protein sequence
3.	Sample.sequencegene.hap2.prot.fa : haplotype 2 protein sequence

## Run Examples

Indiv mode: 

`perl Script/haplotype_translation_final.pl --phased-vcf Example/GJB2_CCDS9290.1.HG0048.vcf --output-dir Run_Example --sample-file Example/HG00448.sample --indiv --synonymous --gene GJB2`

Pop mode: 

`time perl Script/haplotype_translation_final.pl --phased-vcf Example/AFM_CCDS3557.1.example.vcf --output-dir Run_pop_Example --sample-file Example/pop_example.sample --pop --synonymous --gene AFM`

You can add `–-synonymous` argument to view the synonymous mutation on the haplotype.

## Result Examples

To see the examples results go to the Web Application :http://med-laennec.univ-brest.fr/GEMPROT/ 

## Contact Information

Tania Cuppens

Email: tania.cuppens@inserm.fr

## Publication
If you have used Gemprot in a scientific publication, we would appreciate citations to the following paper:

Cuppens, T., Ludwig, T. E., Trouvé, P., & Genin, E. (2018). GEMPROT: visualization of the impact on the protein of the genetic variants found on each haplotype. [Bioinformatics](https://doi.org/10.1093/bioinformatics/bty993)
