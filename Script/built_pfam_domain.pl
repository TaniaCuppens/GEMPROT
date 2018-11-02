#!/usr/bin/perl -w
# Extract dna sequence and translate to prot --- 
# Author: Tania Cuppens 
# Created: 27 september 2017

require Exporter;
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use vars qw(%Cmdline $usage);
use File::Basename;
use File::Glob;
# ':bsd_glob';
use Term::ANSIColor;
use Path::Tiny;
use LWP::UserAgent;
use XML::LibXML;
use Clone qw(clone);
use FindBin qw($Bin);
use lib "$Bin";
use conf::config;
my $cmdAndArgs = join (' ', $0, @ARGV);
my $basename = basename $0;
my $rep_script = dirname $0;
my $mode="";
my $usage="
 
		MANDATORY:
		--gene-file	: official gene symbol with one ccds ; tab delimited file (ex: TP53	CCDS11118.1)
	
    OPTIONNAL:
		-h	: show this message and quit

\n";
# Parse command line
#####################
sub CmdLineParser () {
    Getopt::Long::Configure( "no_ignorecase" );
    #case sensistive
	GetOptions( \%Cmdline,
	"gene-file=s",
	# help
		"h"
    ) || die "$!\n";
    if  (defined $Cmdline {
        h
    }
    ) {
        print colored(['yellow'], $usage);
        exit 0;
    }
}
CmdLineParser();
# main
######
my $gene_file="";
my $n=0;
my @col;
my $chr;
my $gene;
my $gene_start;
my $gene_stop;
my @choose_coding;
my @coding_col;
my $gene_name;
my $gene_line;
my $ccds_num;
my @coding;
my @cod_line;
if (!defined ($Cmdline {
    "gene-file"
}
)) {
    die "no gene file\n";
}
else {
    $gene_file = $Cmdline {
        "gene-file"
    };
    delete $Cmdline {
        "gene-file"
    };
}
open GENE, $gene_file;
my @gene_id = <GENE>;
my %dif_hap=();
foreach $gene_line (@gene_id) {
    chomp $gene_line;
    printf $gene_line."\n";
    my @gene_split = split(/\t/, $gene_line);
    next if (!defined  $gene_split[0]);
    $gene_name=$gene_split[0
    ];
    printf $gene_split[1]."\n";
    if (!defined $gene_split[1]) {
        print "pas defini";
        open(CCDS, "grep -w \"".$gene_name."\" ".$conf::config::CCDS_file." |");
        my @coding = <CCDS>;
        if (@coding == 0
            ){
            print "\n################ \n No coding sequence known for this gene, verify if official gene symbol is correct\n################\n\n";
            exit;
        }
        print "Available coding sequence\n";
        for my $ccds_line (@coding) {
            @col= split(/\t/, $ccds_line);
            my $ccds_id=$col[4];
            open(NM, "grep ".$ccds_id." ".$conf::config::CCDS_NM_file."\n |");
            my @transcript = <NM>;
            print "[".$n."] : ".$ccds_id."  ======>   ";
            for my $nm_line (@transcript) {
                #print $nm_line;
                my @c= split(/\t/, $nm_line);
                my $nm_id=$c[4];
                #print $nm_id;
                if ($nm_id =~ /^NM/) {
                    print $nm_id." ; ";
                }
            }
            print "\n";
            $n++;
        }
        print "Please choose one ccds : ";
        my $inc=<STDIN>;
        chop($inc);
        print $ inc;
        if (looks_like_number($inc) && $inc < $n) {
            @choose_coding=$coding[$inc];
            my @column= split(/\t/, $choose_coding[0]);
            $chr=$column[0
            ];
            $gene_start=$column[7];
            $gene_stop=$column[8];
        }
        else {
            @choose_coding=@coding;
            my @col_umn= split(/\t/, $choose_coding[0]);
            $chr=$col_umn[0
            ];
            $gene_start=$col_umn[7];
            $gene_stop=$col_umn[8];
            print "\n################ \nNot available coding sequence. Using all ccds \n################\n\n";
        }
    }
    else {
        print "defini";
        $ccds_num=$gene_split[1];
        open(CCDS, "grep -w \"".$gene_name."\" ".$conf::config::CCDS_file." | grep -w \"".$ccds_num."\" | ");
        @coding = <CCDS>;
        @cod_line = split(/\t/, $coding[0]);
        $chr=$cod_line[0
        ];
        print Dumper (@coding);
        $gene_start=$cod_line[7];
        $gene_stop=$cod_line[8];
        @choose_coding = @coding;
        if (@coding == 0
            ){
            print "\n################ \n No coding sequence known for this gene, verify if official gene symbol is correct\n################\n\n";
            next;
        }
    }
    foreach my $line (@choose_coding) {
        printf "\ngene ".$line."\n";
        my $sequence;
        my $date_pfam="";
        my $date_release="";
        my @fich1 = split(/\t/, $line);
        $gene=$fich1[2]."_".$fich1[4];
        my $cache=$rep_script."/cache";
        print "\n####Gene name ".$gene."\n";
        if (!-e $cache) {
            mkdir $cache or die "can't create ".$cache."\n";
        }
        my $pfam_store_file = $cache."/pfam_store.txt";
        my $gene_seq_file = $cache."/".$gene.".fa";
        my $pfam="";
        while ($pfam !~ m/^\</) {
            $pfam=`curl --silent -LH 'Expect:' -F output=xml 'http://pfam.xfam.org/family/Piwi/acc'`
            ;
			#print $pfam."\r";
        }
        my $xml_parser = XML::LibXML->new();
        my $dom = $xml_parser->parse_string( $pfam );
        my $root = $dom->documentElement();
        #print $root;
        $date_release = $root->getAttribute( 'release_date' );
        chomp $date_release;
        if (!-e $pfam_store_file) {
            open PSF, ">>".$pfam_store_file;
            print "Existe pas Release_date= ".$date_release ."\n";
            print PSF "Release_date=".$date_release."\n";
            print PSF "Gene\tDomaine\tStart\tEnd\tAccession\n";
        }
        else {
            open(FL, "grep \"Release\" ".$pfam_store_file." |");
            my @firstline =<FL>;
            print "premiere ligne ".Dumper (\@firstline)."\n";
            my @fline=split(/=/, $firstline[0]);
            print Dumper (\@fline);
            $date_pfam=$fline[1];
            chomp $date_pfam;
            if ($date_pfam ne $date_release) {
                print 'Existe mais pas bon Release_date= '.$date_release ."\n";
                system("rm ".$pfam_store_file);
                print "Release_date= ".$date_release ."\n";
                print PSF "Release_date=".$date_release."\n";
                print PSF "Gene\tDomaine\tStart\tEnd\tAccession\n";
            }
        }
        chomp($line);
        my $nb_exon;
        my $nb_aa;
        my $strand=$fich1[6];
        my %h=();
        my $aa_pos;
        my $nt_pos;
        my $gene_start=$fich1[7];
        my $gene_stop=$fich1[8];
        my @exons = split(/,/, $fich1[9]);
        my $nt_nb = 0;
		
        foreach my $exon (@exons) {
            $exon =~ s/[\[\]\s]//g;
            my @pos= split(/-/, $exon);
            my $start=$pos[0]+1;
            my $stop=$pos[1]+1;
            $nb_exon+=1;
            # Keep exon sequence 
			######################
            open(SAMTOOLS, $conf::config::samtools." faidx ".$conf::config::genome_reference." ".$chr.":".$start."-".$stop." | grep -v \">\" | tr -d \"\\n\" |");
            while (my $seqall=<SAMTOOLS>) {
                for my $s (split "", $seqall) {
                    $s=~ tr/atgc/ATGC/;
                    #print $nt_nb."      ".$aa_pos."\n"; 
                    $h {
                        $chr
                    } {
                        $start
                    }
                    =$s;
                    $start++;
                    $nt_nb++;
                }
            }
            close SAMTOOLS;
        }
        if ($fich1[6]eq"+") {
            for my $k( sort {
                $a<=>$b
            }
            keys % {
                $h {
                    $chr
                }
            }
            ) {
                $sequence.=$h {
                    $chr
                } {
                    $k
                };
            }
        }
        elsif ($fich1[6]eq"-") {
            for my $k( sort {
                $a<=>$b
            }
            keys % {
                $h {
                    $chr
                }
            }
            ) {
                my $current_nt1 = $h {
                    $chr
                } {
                    $k
                };
                $current_nt1 =~ tr/ATGCatgc/TACGtacg/;
                $sequence=$current_nt1.$sequence;
            }
        }
        print $sequence;
        my $protein='';
        my $codon;
        for(my $i=0;$i<(length($sequence)-2);
        $i+=3) {
            $codon=substr($sequence,$i,3);
            $protein.=&codon2aa($codon);
        }
        print "The translated protein is :\n$protein\n";
        open (SEQ, ">".$gene_seq_file);
        print SEQ $protein;
        open(GGF, "grep -w \"".$gene."\" ".$pfam_store_file." |");
        my @genedomain =<GGF>;
        if (@genedomain== 0
            ){
            get_url_pfam($gene_seq_file, $pfam_store_file,$gene);
        }
    }
}
sub codon2aa {
    my($codon)=@_;
    $codon=uc $codon;
    my(%g)=(
		'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
		'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L',
		'CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H',
		'CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I',
		'ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
		'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A',
		'GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G',
		'GGG'=>'G','GGT'=>'G'
    );
    if(exists $g {
        $codon
    }
    ) {
        return $g {
            $codon
        };
    }
    else {
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}
# Get pfam result to hash table 
############

sub get_url_pfam {
    my ($sequence_gene_prot,$pfile,$gene) = @_;
    open (PFILE, ">>".$pfile);
    my @url;
    my @result;
    while (@result == 0
        ){
        @url=`curl --silent -LH 'Expect:' -F seq='<$sequence_gene_prot' -F output=xml 'http://pfam.xfam.org/search/sequence'`
        ;
		#print "curl --silent -LH 'Expect:' -F seq='<$sequence_gene_prot' -F output=xml 'http://pfam.xfam.org/search/sequence'\n";
        @result = grep(/result_url/, @url);
        #print Dumper (\@url);
		#print Dumper (\@result);
    }
    $result[0] =~ s/result_url|\>|\<| //g;
    $result[0] =~ s/\/$//;
    my $res_url= $result[0
    ];
    chomp $res_url;
    #$res_url = "'".$res_url."'";
    my $pfam="";
    print "\n==> RUN PFAM\n";
    while ($pfam !~ m/^\</) {
        $pfam=`curl --silent -LH 'Expect:' $res_url`
        ;
        print $pfam."\r";
    }
    print "\n==> END PFAM\n";
    if ($pfam =~m/matches/) {
        #print $pfam;
        my $xml_parser = XML::LibXML->new();
        my $dom = $xml_parser->parse_string($pfam);
        my $root = $dom->documentElement();
        #print Dumper (\$root)."########/n";

        my ( $entry ) = $root->getChildrenByTagName( 'results' );
        my ( $entry1 ) = $entry->getChildrenByTagName( 'matches' );
        my ( $entry2 ) = $entry1->getChildrenByTagName( 'protein' );
        my ( $entry3 ) = $entry2->getChildrenByTagName( 'database' );
        my ( @entry4 ) = $entry3->getChildrenByTagName( 'match' );
        my $ID;
        my $start;
        my $end;
        my $domain;
        my $accession;
        #print Dumper \@entry4;
        my $domain_gene="";
        foreach my $match (@entry4) {
            #print 'domain: ' . $match->getAttribute( 'id' ) . "\n";
            my ( @entry5 ) = $match->getChildrenByTagName( 'location' );
            foreach my $loc (@entry5) {
                $start=$loc->getAttribute( 'ali_start' );
                $end = $loc->getAttribute( 'ali_end' );
                $domain = $match->getAttribute( 'id' );
                $accession = $match->getAttribute( 'accession' );
                $ID++;
                #"Gene\tDomaine\tStart\tEnd\tAccession\n";
                $domain_gene=$domain_gene.$gene."\t".$domain."\t".$start."\t".$end."\t".$accession."\n";
            }
        }
        print PFILE $domain_gene;
    }
}
