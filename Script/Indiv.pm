#!/usr/bin/perl -w
# Population mode ---
# Author: Tania Cuppens
# Created: 23 march 2016

package Indiv;
require Exporter;
use strict;
use warnings;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use LWP::UserAgent;
use XML::LibXML;
use Clone qw(clone);
use experimental 'smartmatch';
sub Indiv {
    my ($arg)=@_;
    my %Cmdline=%$arg;
    # main
	######
    my $in_phased="";
    my $out="";
    my $sample_file="";
    my $gene_file="";
    my $g="";
    my $s="";
    my $domain_file="";
    my $phas="";
    my $phased_vcf="";
    my $syn_opt="";
    my $fa_opt="";
    my $domain_opt="";
    my $pdv="";
    my $rep_script = dirname $0;
    if (!defined ($Cmdline {
        "output-dir"
    }
    )) {
        die "no output directory\n";
    }
    else {
        $out=$Cmdline {
            "output-dir"
        };
        delete $Cmdline {
            "output-dir"
        };
    }
    if (!-e $out) {
        print "Create : $out\n";
        mkdir $out or die "can't create ".$out."\n";
    }
    my $tmp_dir=$out."/tmp/";
    if (!-e $tmp_dir) {
        print "Create : $tmp_dir\n";
        mkdir $tmp_dir or die "can't create ".$tmp_dir."\n";
    }
    my $cache=$rep_script."/cache";
    if (!-e $cache) {
        mkdir $cache or die "can't create ".$cache."\n";
    }
    if (!defined ($Cmdline {
        "phased-vcf"
    }
    )) {
        die "########\nno input file\n########\n";
    }
    else {
        if (-e $Cmdline {
            "phased-vcf"
        }
        && $Cmdline {
            "phased-vcf"
        }
        =~ /\.vcf*/) {
            $in_phased=$Cmdline {
                "phased-vcf"
            };
            delete $Cmdline {
                "phased-vcf"
            };
            if ($in_phased =~ /.gz$/) {
                open(PDV, "zcat ".$in_phased." | head -10 | grep \"fileformat=VCF\" |");
                $pdv = <PDV>;
                $pdv.="";
                if ($pdv eq "") {
                    die "########\nProbably input is not a vcf file\n########\n";
                }
            }
            else {
                open(PDV, "head -10 ".$in_phased." | grep \"fileformat=VCF\" |");
                $pdv = <PDV>;
                $pdv.="";
                if ($pdv eq "") {
                    die "########\nProbably input is not a vcf file\n########\n";
                }
            }
        }
        else {
            die "########\nProbably input is not a vcf file\n########\n";
        }
    }
    my $sample_rec = $tmp_dir."recovery.sample";
    Fonction::get_sample($in_phased,$sample_rec);
    my ($with_chr) = Fonction::get_chr($in_phased);
    if (!defined ($Cmdline {
        "sample-file"
    }
    )) {
        if (!defined ($Cmdline {
            "sample"
        }
        )) {
            printf "#######\nsample recovery in vcf-file\n########\n";
            $sample_file=$sample_rec;
        }
        else {
            $s=$Cmdline {
                "sample"
            };
            #printf $s." sample";
            $sample_file=$tmp_dir.$s.".sample";
            open SAMP, ">".$sample_file;
            print SAMP $s;
            close SAMP;
        }
    }
    else {
        $sample_file = $Cmdline {
            "sample-file"
        };
        #printf $sample_file;
        delete $Cmdline {
            "sample-file"
        };
        if (!-e $sample_file) {
            die "########\nSample-file does not exist\n########\n";
        }
    }
    if (!defined ($Cmdline {
        "gene-file"
    }
    )) {
        if (!defined ($Cmdline {
            "gene"
        }
        )) {
            die "no gene file\n";
        }
        else {
            $g=$Cmdline {
                "gene"
            };
            $gene_file=$tmp_dir.$g.".gene";
            open GENE, ">".$gene_file;
            print GENE $g;
        }
    }
    else {
        $gene_file = $Cmdline {
            "gene-file"
        };
        delete $Cmdline {
            "gene-file"
        };
        if (!-e $gene_file) {
            die "########\nNo such gene-file\n########\n";
        }
    }
    if (!defined ($Cmdline {
        "domain"
    }
    )) {
        $domain_opt = 'no';
    }
    else {
        $domain_file = $Cmdline {
            "domain"
        };
        $domain_opt = 'yes';
        delete $Cmdline {
            "domain"
        };
    }
    if (!defined ($Cmdline {
        "synonymous"
    }
    )) {
        $syn_opt = 'no';
    }
    else {
        $syn_opt = 'yes';
        delete $Cmdline {
            "synonymous"
        };
    }
    if (!defined ($Cmdline {
        "fasta"
    }
    )) {
        $fa_opt = 'no';
    }
    else {
        $fa_opt = 'yes';
        delete $Cmdline {
            "fasta"
        };
    }
    #printf $in_phased." vcf file\n";
	#printf "Fonction::get_sample($in_phased)";

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
    my @non_exist;
    open GENE, $gene_file;
    my @gene_id = <GENE>;
    my %dif_hap=();
    my $index_file=$out."/index.html";
    open INDEX, ">".$index_file;
    print INDEX "<!doctype html>\n<html>";
    print INDEX "<head><title>Choose a transcript</title><link rel=\"stylesheet\" type=\"text/css\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.css\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot16x16.png\" sizes=\"16x16\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot32x32.png\" sizes=\"32x32\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot186x186.png\" sizes=\"186x186\"/><script type=\"text/javascript\" src=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.js\"></script></head>\n";
    print INDEX "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/header.html\" id=\"gemprotheader\" class=\"gemprotheader\"></div>";
    print INDEX "<body>\n<h1>Choose a transcript</h1>\n";
    print INDEX"<div class=\"poptablediv\"><table class=\"poptable\">";
    print INDEX "<thead><tr><th>Gene</th><th>CCDS_ID</th><tr></thead>";
    foreach $gene_line (@gene_id) {
        chomp $gene_line;
        #print $gene_line;
        my @gene_split = split(/\t/, $gene_line);
        next if (!defined  $gene_split[0]);
        $gene_name=$gene_split[0
        ];
        if (!defined $gene_split[1]) {
            open(CCDS, "grep -w \"".$gene_name."\" ".$conf::config::CCDS_file." |");
            my @coding = <CCDS>;
            #print Dumper (@coding);
            if (@coding == 0
                ){
                printf "\n################ \nWARNING : No coding sequence known for $gene_name, verify if this gene symbol officially is correct\n################\n\n";
                exit;
            }
            if (@coding ==1) {
                @choose_coding= @coding;
                @col= split(/\t/, @choose_coding);
                $chr=$col[0
                ];
            }
            else {
                printf "Available coding sequence\n";
                for my $ccds_line (@coding) {
                    @col= split(/\t/, $ccds_line);
                    $chr=$col[0
                    ];
                    my $ccds_id=$col[4];
                    open(NM, "grep ".$ccds_id." ".$conf::config::CCDS_NM_file."\n |");
                    my @transcript = <NM>;
                    printf "[".$n."] : ".$ccds_id."  ======>   ";
                    for my $nm_line (@transcript) {
                        #print $nm_line;
                        my @c= split(/\t/, $nm_line);
                        my $nm_id=$c[4];
                        #print $nm_id;
                        if ($nm_id =~ /^NM/) {
                            printf $nm_id." ; ";
                        }
                    }
                    print "\n";
                    $n++;
                }
                printf "Please choose one ccds : ";
                my $inc;
                alarm(10);
                eval {
                    local $SIG {
                        ALRM
                    }
                    = sub {
                        die
                    };
                    $inc = <STDIN>;
                    alarm 0;
                };
                if ($inc eq "") {
                    printf "time out!\n";
                }
                else {
                    #tludwig is this else needed  ?
                }
                if (looks_like_number($inc) && $inc < $n) {
                    @choose_coding=$coding[$inc];
                }
                else {
                    #printf Dumper (@choose_coding);
                    @choose_coding=@coding;
                    printf "\n################ \nNo available coding sequence. Using all ccds \n################\n\n";
                }
            }
        }
        else {
            $ccds_num=$gene_split[1];
            open(CCDS, "grep -w \"".$gene_name."\" ".$conf::config::CCDS_file." | grep -w \"".$ccds_num."\" | ");
            @coding = <CCDS>;
            #printf Dumper (@coding);
            @choose_coding = @coding;
            if (@coding == 0
                ){
                printf "\n################ \n No coding sequence known for this gene, verify if official gene symbol is correct\n################\n\n";
                next;
            }
        }
        my @cod= split(/\t/, $choose_coding[0]);
        $chr=$cod[0
        ];
        my %clin_hash = &Fonction::clinvartohash($chr, $with_chr);
        # Keep gene transcript sequence
		######################

		#print Dumper (@choose_coding);
        open SAMPLE, $sample_file;
        my @sample_id = <SAMPLE>;
        foreach my $line (@choose_coding) {
            my $nb_sample=0;
			my $nb_hap_tot=0;
            my @fich1 = split(/\t/, $line);
            $chr=$fich1[0
            ];
            if ($with_chr eq "yes") {
                $chr="chr";
            }
            $gene=$fich1[2]."_".$fich1[4];
            $gene_start=$fich1[7];
            $gene_stop=$fich1[8];
            my $gene_dir=$out."/".$gene;
            printf "\n####Gene name ".$gene."\n";
            if (!-e $gene_dir) {
                mkdir $gene_dir or die "can't create ".$gene_dir."\n";
            }
            my $V_dir=$gene_dir."/Visual";
            if (!-e $V_dir) {
                print "Create : $V_dir\n";
                mkdir $V_dir or die "can't create ".$V_dir."\n";
            }
            my $F_dir=$gene_dir."/Fasta";
            if (!-e $F_dir) {
                print "Create : $F_dir\n";
                mkdir $F_dir or die "can't create ".$F_dir."\n";
            }
            print INDEX "<thead><tr><td class=\"gene\">$fich1[2]</td><td><a href=$gene>$fich1[4]</a></td>";
            my $select_phased =$tmp_dir.$gene.".vcf";
            printf "\n==> RUN VCFTOOLS \n";
            if ($in_phased =~ /.gz$/) {
                printf "gzip vcf \n";
                #print $conf::config::vcftools." --gzvcf ".$in_phased." --recode --from-bp ".$gene_start." --to-bp ".$gene_stop." --chr ".$chr." --out ".$select_phased;
                system($conf::config::vcftools." --gzvcf ".$in_phased." --recode --from-bp ".$gene_start." --to-bp ".$gene_stop." --chr ".$chr." --out ".$select_phased." > ".$tmp_dir."/vcftools.out 2> ".$tmp_dir."/vcftools.err");
            }
            else {
                #print $conf::config::vcftools." --vcf ".$in_phased." --recode --from-bp ".$gene_start." --to-bp ".$gene_stop." --chr ".$chr." --out ".$select_phased;
                system($conf::config::vcftools." --vcf ".$in_phased." --recode --from-bp ".$gene_start." --to-bp ".$gene_stop." --chr ".$chr." --out ".$select_phased." > ".$tmp_dir."/vcftools.out 2> ".$tmp_dir."/vcftools.err");
            }
            printf "\n==> END VCFTOOLS \n";
            $phased_vcf=$select_phased.".recode.vcf";
            my $htlm_sum_file=$gene_dir."/index.html";
            open HTMLS, ">".$htlm_sum_file;
            print HTMLS "<!doctype html>\n<html>";
            print HTMLS "<head><title>Summary results </title><link rel=\"stylesheet\" type=\"text/css\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.css\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot16x16.png\" sizes=\"16x16\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot32x32.png\" sizes=\"32x32\"/><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot186x186.png\" sizes=\"186x186\"/><script type=\"text/javascript\" src=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.js\"></script></head>\n";
            print HTMLS "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/header.html\" id=\"gemprotheader\" class=\"gemprotheader\"></div>";
            print HTMLS "<div class=\"goback\"><a href=\"../\">Change Transcript</a></div>";
            print HTMLS "<body>\n<h1>Results summary </h1>\n";
            print HTMLS "<h2> Summary by sample(s)  </h2>\n";
            print HTMLS "<div class=\"freqsummarytablediv\"><table class=\"freqsummarytable\">";
            print HTMLS "<thead><tr><th>Sample</th><th>Haplotype 1</th><th>Haplotype 2</th></tr></thead><tbody>";
            my $gene_mut_file = $tmp_dir."/hap_".$gene.".txt";
            open GMF, ">".$gene_mut_file;
            print GMF "Gene\thaplotype\tnb\tfrequency\n";
            my $bed_file = $tmp_dir."/".$gene.".bed";
            open BED, ">".$bed_file;
            chomp($line);
            my $nb_exon;
            my $nb_aa;
            my $strand=$fich1[6];
            my %h=();
            my %h1=();
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
                print BED $chr."\t".$start."\t".$stop."\n";
                $nb_exon+=1;
                # Keep exon sequence
				######################
                open(SAMTOOLS, $conf::config::samtools." faidx ".$conf::config::genome_reference." ".$fich1[0].":".$start."-".$stop." 2> ".$tmp_dir."/samtools.err | grep -v \">\" | tr -d \"\\n\" |");
                while (my $seqall=<SAMTOOLS>) {
                    for my $s (split "", $seqall) {
                        $s=~ tr/atgc/ATGC/;
                        $aa_pos = int($nt_nb/3)+1;
                        #print $nt_nb."      ".$aa_pos."\n";
                        $h {
                            $chr
                        } {
                            $start
                        } {
                            "nt"
                        }
                        = $s;
                        $h {
                            $chr
                        } {
                            $start
                        } {
                            "amino"
                        }
                        = $aa_pos;
                        $h {
                            $chr
                        } {
                            $start
                        } {
                            "exon"
                        }
                        = $nb_exon;
                        $h {
                            $chr
                        } {
                            $start
                        } {
                            "nt_pos"
                        }
                        = $nt_nb;
                        $h1 {
                            $chr
                        } {
                            $start
                        }
                        = $s;
                        $start++;
                        $nt_nb++;
                    }
                }
                close SAMTOOLS;
            }
            my $sequence_gene_ref=$F_dir."/sequence_".$fich1[4]."_".$gene."_ref_prot.fa";
            my $sequence_name="\">".$fich1[2].".".$chr.":".$fich1[7]."-".$fich1[8]."_".$strand."\"";
            my $dom_file=$tmp_dir."/".$gene.".domain_file.txt";
            my $sequence_protref='';
            my $sequence_protref1='';
            my $sequence_protref2='';
            my $sequence_prot1='';
            my $sequence_prot2='';
            my $sequence_ref='';
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
                    $sequence_ref.=$h {
                        $chr
                    } {
                        $k
                    } {
                        "nt"
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
                    my $current_nt = $h {
                        $chr
                    } {
                        $k
                    } {
                        "nt"
                    };
                    $h {
                        $chr
                    } {
                        $k
                    } {
                        "amino"
                    }
                    =$aa_pos-($h {
                        $chr
                    } {
                        $k
                    } {
                        "amino"
                    }
                    )+1;
                    #print Dumper (\$h{$chr}{$k}{"amino"})."    #####\n";
                    $current_nt =~ tr/ATGCatgc/TACGTACG/;
                    $sequence_ref=$current_nt.$sequence_ref;
                }
            }
            $sequence_protref = Fonction::translate($sequence_ref, $sequence_gene_ref, $sequence_name,$fa_opt);
            if ($sequence_protref eq "" || $sequence_protref eq "*") {
                print "\n################ \n It is non-coding gene. Please check if you taking the data on the same reference genome\n################\n\n";
                next;
            }
            #printf $sequence_gene_ref;
			#printf $sequence_protref;
            my $domain = Fonction::get_url_pfam($sequence_gene_ref, $gene, $cache);
            Fonction::domain_file($dom_file,$domain,$gene);
            my %prot_domain =%$domain;
            my @keys = keys %prot_domain;
            my $nb_domain = @keys;
            my @doma =();
            my %comptage = ();
            foreach my $dom (@keys) {
                push(@doma,$prot_domain {
                    $dom
                } {
                    Domain
                }
                );
            }
            foreach my $mot ( @doma ) {
                $comptage {
                    $mot
                }
                ++;
            }
            #print Dumper (\@doma);
            foreach my $line_sample (@sample_id) {
                chomp($line_sample);
                next if ($line_sample=~/^$/);
                my @sep_sample= split(/\t/, $line_sample);
                my $sample= $sep_sample[0
                ];
                open(PVCF, "grep -w \"".$sample."\" ".$phased_vcf." |");
                my @exist_sample = <PVCF>;
                if (@exist_sample == 0
                    ){
                    printf "sample $sample doesn't exist in the vcf\n";
                    push(@non_exist,$sample);
                    next;
                }
                my %hashap1 = % {
                    clone(\%h1)
                };
                #print Dumper (\%h1);
                my %hashap2 = % {
                    clone(\%h1)
                };
                my %hasref1 = % {
                    clone(\%h1)
                };
                my %hasref2 = % {
                    clone(\%h1)
                };
                #print Dumper (\%h2);
                $nb_aa = $aa_pos-1;
                my $htlm_file=$V_dir."/".$sample."_".$gene."_result.html";
                open HTML, ">".$htlm_file;
                print HTML "<!doctype html>\n<html>\n<head>\n<title>$sample result - $gene Gene</title><link rel=\"stylesheet\" type=\"text/css\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.css\"><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot16x16.png\" sizes=\"16x16\"><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot32x32.png\" sizes=\"32x32\"><link rel=\"icon\" type=\"image/png\" href=\"http://lysine.univ-brest.fr/GEMPROT/gemprot186x186.png\" sizes=\"186x186\"><script type=\"text/javascript\" src=\"http://lysine.univ-brest.fr/GEMPROT/gemprot.js\"></script>\n</head>\n";
                print HTML "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/header.html\" id=\"gemprotheader\" class=\"gemprotheader\"></div>";
                print HTML "<div class=\"goback\"><a href=\"../\">Results summary</a></div>";
                print HTML "<body>\n<h1>$sample result - $gene Gene</h1>\n";
                print HTML"<p class=\"geneinfo\"> $gene is a protein with $nb_aa amino acids and $nb_exon exon(s) <br> The protein has $nb_domain domain(s) :<br>";
                my @dom_list;
                #my $dom="";

                foreach my $domb (sort keys (%prot_domain)) {
                    my $name = $prot_domain {
                        $domb
                    } {
                        "Domain"
                    };
                    if (grep {
                        $_ eq $name
                    }
                    @dom_list) {
                        #tluldwig reverted grep ?
                    }
                    else {
                        print HTML  $comptage {
                            $prot_domain {
                                $domb
                            } {
                                "Domain"
                            }
                        }
                        ." <a href=\"http://pfam.xfam.org/family/".$prot_domain {
                            $domb
                        } {
                            "Accession"
                        }
                        ."\">".$prot_domain {
                            $domb
                        } {
                            "Domain"
                        }
                        ."</a> domain(s) <br>";
                        push(@dom_list,$prot_domain {
                            $domb
                        } {
                            "Domain"
                        }
                        );
                    }
                }
                #print HTML"<p class=\"sampleresult\"> <h2>Results for $sample.</h2>\n<br>";
                chomp($sample);
                $nb_sample+=1;
                #print "\n\n==> $sample\n\n";
                print "Sample ".$nb_sample."   \r";
                # Phased and add sample mutation
				######################

                my $sequence1='';
                my $sequence2='';
                my $sequenceref1='';
                my $sequenceref2='';
                my %hap1 = ();
                my %hap2 = ();
                my %ref1 = ();
                my %ref2 = ();
                my %mut_hap = ();
                my %mut_hap_aa = ();
                my ($haplotype1, $haplotype2, $referencehap1, $referencehap2, $mut_hash, $mut_hash_pos,$nb1, $nb2) = Fonction::parsevcf2($phased_vcf,\%h,\%hashap1,\%hashap2,\%hasref1,\%hasref2,$sample,\%clin_hash,$gene);
                %hap1=%$haplotype1;
                %hap2=%$haplotype2;
                %ref1=%$referencehap1;
                %ref2=%$referencehap2;
                %mut_hap=%$mut_hash;
                %mut_hap_aa=%$mut_hash_pos;
                #print Dumper (\%mut_hap);
				#print Dumper (\%mut_hap_aa);
				
                if ($fich1[6]eq"+") {
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $hap1 {
                            $chr
                        }
                    }
                    ) {
                        $sequence1.=$hap1 {
                            $chr
                        } {
                            $k
                        };
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $hap2 {
                            $chr
                        }
                    }
                    ) {
                        $sequence2.=$hap2 {
                            $chr
                        } {
                            $k
                        };
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $ref1 {
                            $chr
                        }
                    }
                    ) {
                        $sequenceref1.=$ref1 {
                            $chr
                        } {
                            $k
                        };
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $ref2 {
                            $chr
                        }
                    }
                    ) {
                        $sequenceref2.=$ref2 {
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
                        $hap1 {
                            $chr
                        }
                    }
                    ) {
                        my $current_nt1 = $hap1 {
                            $chr
                        } {
                            $k
                        };
                        $current_nt1 =~ tr/ATGCatgc/TACGtacg/;
                        $sequence1=$current_nt1.$sequence1;
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $hap2 {
                            $chr
                        }
                    }
                    ) {
                        my $current_nt2 = $hap2 {
                            $chr
                        } {
                            $k
                        };
                        $current_nt2 =~ tr/ATGCatgc/TACGtacg/;
                        $sequence2=$current_nt2.$sequence2;
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $ref1 {
                            $chr
                        }
                    }
                    ) {
                        my $current_nt2 = $ref1 {
                            $chr
                        } {
                            $k
                        };
                        $current_nt2 =~ tr/ATGCatgc/TACGtacg/;
                        $sequenceref1=$current_nt2.$sequenceref1;
                    }
                    for my $k( sort {
                        $a<=>$b
                    }
                    keys % {
                        $ref2 {
                            $chr
                        }
                    }
                    ) {
                        my $current_nt2 = $ref2 {
                            $chr
                        } {
                            $k
                        };
                        $current_nt2 =~ tr/ATGCatgc/TACGtacg/;
                        $sequenceref2=$current_nt2.$sequenceref2;
                    }
                }
                my $sequence_gene_prot1=$F_dir."/".$sample.".sequence_".$gene."_hap1_prot.fa";
                my $sequence_gene_prot2=$F_dir."/".$sample.".sequence_".$gene."_hap2_prot.fa";
                my $sequence_gene_refbis=$F_dir."/".$sample.".sequence_".$gene."_refv_prot.fa";
                # Translation
				
				######################
                $sequence_prot1 = Fonction::translate($sequence1, $sequence_gene_prot1, $sequence_name,$fa_opt);
                #print " ###########sequence 1 \n".$sequence1." \n";

                $sequence_prot2 = Fonction::translate($sequence2, $sequence_gene_prot2, $sequence_name,$fa_opt);
                #print " ###########sequence 1 \n".$sequence1." \n";
                $sequence_protref1 = Fonction::translate($sequenceref1, $sequence_gene_refbis, $sequence_name,$fa_opt);
                #print " ###########sequence 1 \n".$sequenceref1." \n";
                $sequence_protref2= Fonction::translate($sequenceref2, $sequence_gene_refbis, $sequence_name,$fa_opt);
                #print " ###########sequence 1 \n".$sequenceref2." \n";
				#my ($sequence_ALref1,$sequence_ALprot1)=Fonction::Sequence_comparaison($sequence_protref1, $sequence_prot1);
				#my ($sequence_ALref2,$sequence_ALprot2)=Fonction::Sequence_comparaison($sequence_protref2, $sequence_prot2);

                system("rm ".$sequence_gene_refbis);
                my $mut_1;
                my $mut_2;
                my %hash_mut_type1;
                my %hash_mut_type2;
                my ($mut_info1, $mut_1b,  $mut_type1)=Fonction::prot_modif_type($sequence_protref1,$sequence_prot1,\%mut_hap_aa,$syn_opt,$chr,'hap1');
                my ($mut_info2, $mut_2b, $mut_type2)=Fonction::prot_modif_type($sequence_protref2,$sequence_prot2,\%mut_hap_aa,$syn_opt,$chr,'hap2');
                %hash_mut_type1=%$mut_type1;
                %hash_mut_type2=%$mut_type2;
                #printf $mut_1b." mut_1b\n";
				#printf $mut_2b." mut_2b\n";
				#printf $sample." sample\n";

                if (defined $dif_hap {
                    $gene
                } {
                    $mut_1b
                }
                ) {
                    $dif_hap {
                        $gene
                    } {
                        $mut_1b
                    } {
                        "nb"
                    }
                    +=1;
                    $dif_hap {
                        $gene
                    } {
                        $mut_1b
                    } {
                        "Sample"
                    }
                    .= ",".$sample;
                    $nb_hap_tot+=1;
                }
                else {
                    $dif_hap {
                        $gene
                    } {
                        $mut_1b
                    } {
                        "nb"
                    }
                    =1;
                    $dif_hap {
                        $gene
                    } {
                        $mut_1b
                    } {
                        "Sample"
                    }
                    = $sample;
                    $nb_hap_tot+=1;
                }
                if (defined $dif_hap {
                    $gene
                } {
                    $mut_2b
                }
                ) {
                    $dif_hap {
                        $gene
                    } {
                        $mut_2b
                    } {
                        "nb"
                    }
                    +=1;
                    $dif_hap {
                        $gene
                    } {
                        $mut_2b
                    } {
                        "Sample"
                    }
                    .= ",".$sample;
                    $nb_hap_tot+=1;
                }
                else {
                    $dif_hap {
                        $gene
                    } {
                        $mut_2b
                    } {
                        "nb"
                    }
                    =1;
                    $dif_hap {
                        $gene
                    } {
                        $mut_2b
                    } {
                        "Sample"
                    }
                    = $sample;
                    $nb_hap_tot+=1;
                }
                print HTMLS "<tr><td><a href=Visual/".$sample."_".$gene."_result.html>".$sample."</a></td><td class=\"haplotype\">$mut_1b</td><td class=\"haplotype\">$mut_2b</td>";
                #print Dumper (\%dif_hap);
                my $sample_mut_file = $tmp_dir.$gene.".mutation_by_hap_".$sample.".txt";
                open SMF, ">".$sample_mut_file;
                print SMF "Hap\tGene\tPosition\tRef\tAlt\tType\tStyle\tColor\n";
                if ($mut_info1 ne "no_mutations") {
                    for my $m (split ",", $mut_info1) {
                        my @lmut = split(/:/, $m);
                        #print $mut_1.",hap1\t".$gene."\t".$lmut[1]."\t".$lmut[0]."\t".$lmut[2]."\n";
                        print SMF $mut_info1.",hap1\t".$gene."\t".$lmut[1]."\t".$lmut[0]."\t".$lmut[2]."\t".$lmut[3]."\t".$lmut[4]."\t".$lmut[5]."\n";
                    }
                }
                if ($mut_info2 ne "no_mutations") {
                    for my $m2 (split ",", $mut_info2) {
                        my @lmut2 = split(/:/, $m2);
                        #print $mut_2.",hap2\t".$gene."\t".$lmut2[1]."\t".$lmut2[0]."\t".$lmut2[2]."\n";
                        print SMF $mut_info2.",hap2\t".$gene."\t".$lmut2[1]."\t".$lmut2[0]."\t".$lmut2[2]."\t".$lmut2[3]."\t".$lmut2[4]."\t".$lmut2[5]."\n";
                    }
                }
                system("Rscript ".$rep_script."/Plot_new_annotation.R -m ".$sample_mut_file." -d ".$dom_file." -l ".$aa_pos." -p ".$V_dir."/".$gene."_".$sample." -s yes > /dev/null 2>&1");
                #print "Rscript Plot_new_annotation.R -m ".$sample_mut_file." -d ".$dom_file." -l ".$aa_pos." -p ".$V_dir."/".$gene."_".$sample." -s yes\n";
                my $pdf= $gene."_".$sample."_protein_plot.pdf";
                Fonction::html_file($chr, $aa_pos, $htlm_file, $pdf, $nb1, $nb2, $syn_opt, \%mut_hap, \%hash_mut_type1, \%hash_mut_type2, \%clin_hash);
            }
            print HTMLS "</tbody></table></div>";
            #printf Dumper (\%dif_hap);
			# hap by pop file
			######################

            my $hm;
            foreach $hm (keys % {
                $dif_hap {
                    $gene
                }
            }
            ) {
                #print $gene."\t".$hm."\t";
                print GMF $gene."\t".$hm."\t".$dif_hap {
                    $gene
                } {
                    $hm
                }
                ."\t".$dif_hap {
                    $gene
                } {
                    $hm
                } {
                    "nb"
                }
                /$nb_hap_tot."\t".$dif_hap {
                    $gene
                } {
                    $hm
                } {
                    "Sample"
                }
                ."\n";
                #print HTMLS $gene."\t".$hm."\t".$dif_hap{$gene}{$hm}."\t".$dif_hap{$gene}{$hm}/$nb_hap_tot."\n";
            }
            my $nb_indiv= $nb_hap_tot/2;
            #print HTMLS "</tr>";

            print HTMLS "<h2> Summary by haplotype(s)  </h2>\n";
            print HTMLS "<div class=\"freqsummarytablediv\"><table class=\"freqsummarytable\">";
            print HTMLS "<thead><tr><th>Gene</th><th>Haplotype</th><th>No. of chromosomes</th><th>Frequency</th><th>Individuals</th></tr></thead><tbody>";
            my $hp;
            foreach $hp (keys % {
                $dif_hap {
                    $gene
                }
            }
            ) {
                #print $gene."\t".$hm."\t";
                my $freq = sprintf("%.3g", $dif_hap {
                    $gene
                } {
                    $hp
                } {
                    "nb"
                }
                /$nb_hap_tot);
                print HTMLS "<tr><td class=\"gene\">$gene</td><td class=\"haplotype\">$hp</td><td class=\"nchr\">$dif_hap{$gene}{$hp}{\"nb\"}</td><td class=\"freq\">$freq</td><td>";
                my @sample_split = split(/,/, $dif_hap {
                    $gene
                } {
                    $hp
                } {
                    "Sample"
                }
                );
                foreach my $sp (@sample_split) {
                    print HTMLS "<a href=Visual/".$sp."_".$gene."_result.html>".$sp."</a> ";
                    #<div class=\"box\"><iframe src=\"".$gene."_".$sp."_protein_plot.pdf\" width = \"500px\" height = \"300px\"></iframe></div>";
                }
                print  HTMLS "</td></tr>";
            }
            print HTMLS "</tbody></table></div>";
            print HTMLS "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/footer.html\" id=\"gemprotfooter\" class=\"gemprotfooter\"></div>";
            print HTMLS "<script>includeHTML();</script></body>";
        }
    }
    print INDEX "</tbody></table></div>";
    print INDEX "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/footer.html\" id=\"gemprotfooter\" class=\"gemprotfooter\"></div>";
    print INDEX "<script>includeHTML();</script></body></html>";
    close GMF;
    close CCDS;
    close GENE;
}
1;
