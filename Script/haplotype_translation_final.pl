#!/usr/bin/perl -w
# Extract dna sequence and translate to prot --- 
# Author: Tania Cuppens 
# Created: 23 march 2016

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use vars qw(%Cmdline $usage);
use File::Basename;
use File::Glob;# ':bsd_glob';
use Term::ANSIColor;

use FindBin qw($Bin);
use lib "$Bin";
use conf::config;
use Fonction_final;
use Pop;
use Indiv;

my $cmdAndArgs = join (' ', $0, @ARGV);
my $basename = basename $0;
my $rep_script = dirname $0;
my $mode="";
my $usage="
 
2 modes :
	-indiv
		MANDATORY:
		--phased-vcf	: phased vcf
		--output-dir	: specify output directory
		--sample | --sample-file	: sample or sample file (one sample name by line)
		--gene | --gene-file	: official gene symbol or official gene symbol with one ccds in tab delimited file (ex: TP53	CCDS11118.1) 
	-pop 
		MANDATORY:
		--phased-vcf	: phased vcf
		--output-dir	: specify output directory
		--sample-file	: sample name with his location ; tab delimited file (ex: HG00096	EUR)
		--gene | --gene-file	: official gene symbol or official gene symbol with one ccds in tab delimited file (ex: TP53	CCDS11118.1)
		--loc-file	: location file
\n";
my $opt="
    OPTIONAL:
            --fasta : if you want protein fasta file for each haplotype and reference
		--synonymous	: if you want to see synonymous SNP
            --domain : if you know protein domain 
		-h	: show this message and quit

\n";

# Parse command line
#####################
sub CmdLineParser () {
    Getopt::Long::Configure( "no_ignorecase" ); #case sensistive
    GetOptions( \%Cmdline,
	"indiv",
	"pop",
	"phased-vcf=s",
        "output-dir=s",
        "sample-file=s",
        "sample=s",
        "gene-file=s",
        "gene=s",
	"loc-file=s",
	"synonymous",
 "fasta",
 "domain",
        # help
        "h"
    ) || die "$!\n";
	
    if  (defined $Cmdline{h}) {
        printf colored(['yellow'], $usage);
        printf colored(['grey15'], $opt);
        $Cmdline{"output-dir"} = $rep_script."/cache/";
       exit 0;
    }
	
	if (defined $Cmdline{indiv}){
		$mode = "indiv";
		print "Entering in --".$mode." GEMPROT mode\n";
	}
	
	elsif (defined $Cmdline{pop}){
		$mode = "pop";
		print "Entering in --".$mode." GEMPROT mode\n";
		
	}
	else {
		$mode = "indiv";
		print "Entering in --".$mode." GEMPROT mode\n";;
	
	}
}

CmdLineParser();



if ($mode eq "pop"){
	Pop::Pop(\%Cmdline);
	}
elsif ($mode eq "indiv"){
	Indiv::Indiv(\%Cmdline);
	}
else {
	exit;
}

system("touch ".$Cmdline{"output-dir"}."/tmp/OK");

END{system("touch ".$Cmdline{"output-dir"}."/tmp/OK")};
printf "\nDONE \n";
