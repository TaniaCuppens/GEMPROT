# GEMPROT
Visualization of the impact on the protein of genetic variants found on each haplotype

## Introduction
GEMPROT (Genetic Mutation to Protein Translation) is a two mode tool using phased mutations present within a gene to restore its haplotypes. It produces the protein sequences translated from each chromosome to identify the pattern of mutations in the gene. The first mode outputs both sequences for each individuals and provides a frequency summary of all haplotypes. The second mode is intended for larger populations and shows difference of genes haplotypes repartition.

## Table of content

-	[Dependencies](#Dependencies)
-	[Configuration](#Configuration)
-	[Running GEMPROT](#Running-GEMPROT)
-	[GEMPROT Input](#GEMPROT-Input)
-	[GEMPROT Output](#GEMPROT-Output)
-	[RunExamples](#RunExamples)
-	[Contact Information](#Contact-Information)
-	[License Agreement](#License-Agreement)

## Dependencies

GEMPROT is written in Perl v5.22.1 and requires the following packages before running:
-	[Vcftools](https://sourceforge.net/projects/vcftools/files/) 
-	[Samtools](https://sourceforge.net/projects/samtools/files/samtools/0.1.19/)  
-	R with "optparse" package
`install.packages("optparse")`

GEMPROT run with a reference genome (fasta), with two CCDS files (Consensus Coding Sequence) and with clinvar data file. These files already exist in the download directory but if you want to use others files please change their path in the configuration file. 
You can download in the FTP sites (be careful about taking the data on the same reference genome): 
- ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/  
- ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/  

GEMPROT needs an internet access to import domain information via [Pfam](https://pfam.xfam.org/) (Protein families database).

## Configuration


GEMPROT needs some path information.  

conf/conf.pm : 

Reference
############

$genome_reference= "/FILE_PATH/*.fasta";

$CCDS_file="/FILE_PATH /*.txt";

$CCDS_NM_file="/FILE_PATH /*.txt";

$clin_var="/FILE_PATH /clinvar_*.vcf.gz";

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
        -pop
    MANDATORY:
                --phased-vcf    : phased vcf
                --output-dir    : specify output directory
                --sample-file   : sample name with his location ; tab delimited file (ex: HG00096       EUR)
                --gene | --gene-file    : official gene symbol or official gene symbol with one ccds in tab delimited file (ex: TP53    CCDS11118.1)
                --loc-file      : location file


    OPTIONAL:
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
With –fasta option:
Genename/Fasta/
3.	Sample.sequencegene.ref.prot.fa : reference protein sequence
4.	Sample.sequencegene.hap1.prot.fa : haplotype 1 protein sequence
5.	Sample.sequencegene.hap2.prot.fa : haplotype 2 protein sequence

## Run Examples

Indiv mode: 

`perl Script/haplotype_translation_final.pl --phased-vcf Example/example.vcf --output-dir Results/Run_Example --sample-file Example/sample.txt --gene-file Example/chr1_gene_CCDS.gene –indiv`

Pop mode: 

`perl Script/haplotype_translation_final.pl --phased-vcf Example/example.vcf --output-dir Results/Run_Example --sample-file Example/sample.txt --loc-file Example/loc.txt --gene-file Example/chr1_gene_CCDS.gene --pop`

You can add `–-synonymous` argument to view the synonymous mutation on the haplotype.

## Result Examples

As mentioned 

## Contact Information

Tania Cuppens

Email: tania.cuppens@inserm.fr

## License Agreement

