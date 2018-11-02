#!/usr/bin/perl -w
# Extract dna sequence and translate to prot ---
# Author: Tania Cuppens
# Created: 23 march 2016
package Fonction;
require Exporter;
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use LWP::UserAgent;
use XML::LibXML;
use Clone qw(clone);
use Scalar::Util qw(looks_like_number);
sub checkfileconfig {
    my $refg="";
    my $refi="";
    my $ccdsf="";
    my $ccdsnm="";
    my $clinf="";
    my $samb="";
    my $vcft="";
    open(REFG, "ls -l ".$conf::config::genome_reference."*.fai 2> /dev/null | ");
    $refg.=<REFG>;
    my $ref_dir = dirname $conf::config::genome_reference;
    open(REFI, "ls -l ".$ref_dir."/*.fai 2> /dev/null | ");
    $refi.=<REFI>;
    if ($refg eq "" || $refi eq "") {
        die "\nERROR : reference genome is not found check configuration file, conf.pm \n";
    }
    open(CCDSF, "ls -l ".$conf::config::CCDS_file." 2> /dev/null | ");
    $ccdsf.=<CCDSF>;
    if ($ccdsf eq "") {
        die "\nERROR : ccds file is not found check configuration file, conf.pm \n";
    }
    open(CCDSNM, "ls -l ".$conf::config::CCDS_NM_file." 2> /dev/null | ");
    $ccdsnm.=<CCDSNM>;
    if ($ccdsnm eq "") {
        die "\nERROR : ccds NM file is not found check configuration file, conf.pm \n";
    }
    open(CLINF, "ls -l ".$conf::config::clin_var." 2> /dev/null | ");
    $clinf.=<CLINF>;
    if ($clinf eq "") {
        die "\nERROR : clinvar file is not found check configuration file, conf.pm \n";
    }
    open(SAMTOOLSB, "ls -l ".$conf::config::samtools." 2> /dev/null | ");
    $samb.=<SAMTOOLSB>;
    if ($samb eq "") {
        die "\nERROR : samtools is not found check configuration file, conf.pm \n";
    }
    open(VCFTOOLS, "ls -l ".$conf::config::vcftools." 2> /dev/null | ");
    $vcft.=<VCFTOOLS>;
    if ($vcft eq "") {
        die "\nERROR : vcftools is not found check configuration file, conf.pm \n";
    }
}
sub clinvartohash {
    my ($chr, $with_chr)=@_;
    my @CLNORIGIN;
    my @CLNSIG;
    my @VC;
    my @CLNDBN;
    my @ALT;
    my @IDRS;
    my @info;
    my ($CHROM,	$POS, $ID, $REF, $ALTS, $QUAL, $FILTER, $INFO);
    my $nb_clin_alt;
    #print $chr."\n\n";
    my @clin_chr=`zgrep -w "^$chr" $conf::config::clin_var`
    ;
	#open CLI, "/home/tcuppens/working_directory/Shapeit/update/clintest.vcf";
	#my @clin_chr=<CLI>;
	#print Dumper (\@clin_chr);
    my %table;
    my $nb_alt=0;
    foreach my $l (@clin_chr) {
        #print $l;
        my ($CHROM,	$POS, $ID, $REF, $ALTS, $QUAL, $FILTER, $INFO) = split(/\t/, $l);
        if ($with_chr eq "yes") {
            $CHROM = "chr".$CHROM;
        }
        @ALT=split(/,/, $ALTS);
        @info = split(/;
        /, $INFO);
        if (defined $table {
            $CHROM
        } {
            $POS
        }
        ) {
            #printf $l."\n";
            $nb_alt+=1;
            $table {
                $CHROM
            } {
                $POS
            } {
                "HGVS".$nb_alt
            }
            = {
                ID => $ID,
				REF  => $REF,
				QUAL  => $QUAL,
				FILTER => $FILTER,
				ALT => $ALTS
            };
            foreach my $id_info (@info) {
                if ($id_info=~ /^CLNHGVS/) {
                    my @CLNHGVS = split(/=/, $id_info);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "NC"
                    }
                    =$CLNHGVS[1];
                }
                if ($id_info=~ /^RS=/) {
                    # CLNORIGIN CLNSIG VC CLNDBN
                    @IDRS = split(/=/, $id_info);
                    #print Dumper (\@CLNORIGIN);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "RS"
                    }
                    = "rs".$IDRS[1];
                }
                elsif ($id_info=~ /^CLNORIGIN/) {
                    # CLNORIGIN CLNSIG VC CLNDBN
                    @CLNORIGIN = split(/=/, $id_info);
                    #print Dumper (\@CLNORIGIN);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "CLNORIGIN"
                    }
                    = $CLNORIGIN[1];
                }
                elsif ($id_info=~ /^CLNSIG/) {
                    # CLNORIGIN CLNSIG VC CLNDBN
                    @CLNSIG = split(/=/, $id_info);
                    #print Dumper (\@CLNSIG);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "CLNSIG"
                    }
                    = $CLNSIG[1];
                }
                elsif ($id_info=~ /^CLNDBN/ || $id_info=~ /^CLNDN/) {
                    # CLNORIGIN CLNSIG VC CLNDBN
                    @CLNDBN = split(/=/, $id_info);
                    #print Dumper (\@CLNSIG);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "CLNDBN"
                    }
                    = $CLNDBN[1];
                }
                elsif ($id_info=~ /^VC/) {
                    # CLNORIGIN CLNSIG VC CLNDBN
                    @VC = split(/=/, $id_info);
                    #print Dumper (\@CLNSIG);
                    $table {
                        $CHROM
                    } {
                        $POS
                    } {
                        "HGVS".$nb_alt
                    } {
                        "VC"
                    }
                    = $VC[1];
                }
            }
        }
        else {
            $nb_alt=0;
			my $nb_clin_alt=0;
			#print Dumper (\@info);
            foreach my $id_info (@info) {
                #print $id_info;
                if ($id_info=~ /^CLNHGVS/) {
                    my @CLNHGVS = split(/=/, $id_info);
                    my @CLNHGVS1=split(/\,/, $CLNHGVS[1]);
                    #print Dumper (\@CLNHGVS1);
                    foreach my $id_clnhgvs(@CLNHGVS1) {
                        #print $id_clnhgvs;
						#$table{$CHROM}{$POS}{"HGVS".$nb_clin_alt}{"NC"} =1;#$CLNHGVS1[$nb_clin_alt];
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        }
                        = {
                            NC=> $CLNHGVS1[$nb_clin_alt],
							ID => $ID,
							REF  => $REF,
							QUAL  => $QUAL,
							FILTER => $FILTER,
							ALT => $ALT[$nb_clin_alt]
                        };
                    }
                }
            }
            for (my $i=0;$i<@ALT;
            $i++) {
                foreach my $id_info (@info) {
                    if ($id_info=~ /^RS=/) {
                        # CLNORIGIN CLNSIG VC CLNDBN
                        @IDRS = split(/=/, $id_info);
                        my @IDRS1=split(/\,/, $IDRS[1]);
                        #print Dumper (\@CLNORIGIN);
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        } {
                            "RS"
                        }
                        = "rs".$IDRS1[$nb_clin_alt];
                    }
                    elsif ($id_info=~ /^CLNORIGIN/) {
                        # CLNORIGIN CLNSIG VC CLNDBN
                        @CLNORIGIN = split(/=/, $id_info);
                        my @CLNORIGIN1=split(/\,/, $CLNORIGIN[1]);
                        #print Dumper (\@CLNORIGIN);
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        } {
                            "CLNORIGIN"
                        }
                        = $CLNORIGIN1[$nb_clin_alt];
                    }
                    elsif ($id_info=~ /^CLNSIG/) {
                        # CLNORIGIN CLNSIG VC CLNDBN
                        @CLNSIG = split(/=/, $id_info);
                        my @CLNSIG1=split(/\,/, $CLNSIG[1]);
                        #print Dumper (\@CLNSIG);
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        } {
                            "CLNSIG"
                        }
                        = $CLNSIG1[$nb_clin_alt];
                    }
                    elsif ($id_info=~ /^CLNDBN/ || $id_info=~ /^CLNDN/) {
                        # CLNORIGIN CLNSIG VC CLNDBN
                        @CLNDBN = split(/=/, $id_info);
                        my @CLNDBN1=split(/\,/, $CLNDBN[1]);
                        #print Dumper (\@CLNSIG);
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        } {
                            "CLNDBN"
                        }
                        = $CLNDBN1[$nb_clin_alt];
                    }
                    elsif ($id_info=~ /^VC/) {
                        # CLNORIGIN CLNSIG VC CLNDBN
                        @VC = split(/=/, $id_info);
                        my @VC1=split(/\,/, $VC[1]);
                        #print Dumper (\@CLNSIG);
                        $table {
                            $CHROM
                        } {
                            $POS
                        } {
                            "HGVS".$nb_clin_alt
                        } {
                            "VC"
                        }
                        = $VC1[$nb_clin_alt];
                    }
                }
                $nb_clin_alt+=1;
            }
            #print Dumper (\%table);
        }
    }
    #print Dumper (\%table);
    return %table;
}
# Translation
############
sub translate {
    my($seq,$out,$seq_name,$opt)=@_;
    if ($out =~ "ref") {
        $opt="yes";
    }
    $seq_name=~ s/"//g;
	$seq =~ s/\s//g;
	my $protein='';
	my $codon;
	for(my $i=0;$i<(length($seq)-2);$i+=3){
		my $del;
		my $aa='aa';
		my $pos=$i;
		my $x=0;
		$codon=substr($seq,$i,3);
		#print $codon."\n";
		if ($codon=~m/-/){
			#print "il y a une deletion $codon\n";
			while ($codon =~ /-/g && length($aa)<3 && $pos<(length($seq)-2)) {
				$x++;
				$del++;
				$codon=substr($seq,$i,3+$x);
				$pos = $i+3+$x;
				$aa=$codon;
				$aa =~ s/-//g;
				if ($aa eq "" && $pos>(length($seq)-2)){
					$protein.="-";
					#print "la deletion est de 3 nt\n";
				}
			}
			if ($del/3>=1){
				#print "c'est une deletion de + de 3 nt\n";
				my $j;
				$protein.=&codon2aa($aa);
				for($j=1;$j<=$del/3;$j++){
					$protein.="-";
				}
				$i=$i+$x;
				$del=$del-(($j-1)*3);
			}
			elsif ($pos>(length($seq)-2)){
				my $prot=$protein;
				$prot=~ s/-//g;
				if ($opt eq "yes"){
					open PROTFILE, ">".$out;
					print PROTFILE $seq_name."\n".$prot;
				}
				return $protein."\n";
			}
			else {
				$protein.=&codon2aa($aa);
				$i=$i+$x;
			}
		}
		else {
			$protein.=&codon2aa($codon);
		}
		if ($protein =~ /\*/g){
			my $prot=$protein;
			$prot=~ s/-//g;
			if ($opt eq "yes"){
				open PROTFILE, ">".$out;
				print PROTFILE $seq_name."\n".$prot;
			}
			return $protein."\n";
		}
	}
	my $prot=$protein;
	$prot=~ s/-//g;
	if ($opt eq "yes"){
		open PROTFILE, ">".$out;
		print PROTFILE $seq_name."\n".$prot;
	}
	close PROTFILE; #tludwig close was after return, so it was never executed
	return $protein."\n";
}
sub AlleleOrigin{
	my($id_origin)=@_;
	$id_origin=uc $id_origin;
	my(%o)=(
		'0' => 'unknown',
		'1' => 'germline',
		'2' => 'somatic',
		'3' => '3?',
		'4' => 'inherited',
		'8' => 'paternal',
		'9' => '9?',
		'16' => 'maternal',
		'32' => 'de-novo',
		'64' => 'biparental',
		'128' => 'uniparental',
		'256' => 'not-tested',
		'512' => 'tested-inconclusive',
		'1073741824' => 'other'
    );
    if(exists $o {
        $id_origin
    }
    ) {
        return $o {
            $id_origin
        };
    }
    else {
        return "not in info origin";
    }
}
sub VariantClinicalSignificance {
    my($id_sign)=@_;
    $id_sign=uc $id_sign;
    my(%c)=(
		'0' => 'Uncertain significance',
		'1' => 'not provided',
		'2' => 'Benign',
		'3' => 'Likely benign',
		'4' => 'Likely pathogenic',
		'5' => 'Pathogenic',
		'6' => 'drug response',
		'7' => 'histocompatibility',
		'255' => 'other'
    );
    if (looks_like_number($id_sign) && exists $c {
        $id_sign
    }
    ) {
        return $c {
            $id_sign
        };
    }
    else {
        return $id_sign;
    }
}
sub codon2aa {
    my($codon)=@_;
    $codon=uc $codon;
    my(%g)=(
		'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
		'TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W',
		'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
		'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
		'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
		'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
		'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
		'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G'
    );
    if(exists $g {
        $codon
    }
    ) {
        return $g {
            $codon
        };
    }
    elsif($codon=~m/,/) {
        print STDERR "\n########\n verify your input, only work on bi-allelic SNP and indel file !! \n########\n";
        exit;
    }
    else {
        print STDERR "\n########\n Bad codon \"$codon\"!!\n########\n";
        exit;
    }
}
sub aa2properties {
    my($aa)=@_;
    #chomp $aa;
    $aa=uc $aa;
    my(%p)=('G'=>'smallest amino acid, aliphatic, no charge, not hydrophilic (0.67)', 'A'=>'no charge, aliphatic, hydrophobic (1.0)' , 'V'=>'no charge, aliphatic, hydrophobic (2.3) ', 'L'=>'no charge, aliphatic, hydrophobic (2.2), isomer of isoleucine (I)', 'I'=>'no charge, aliphatic, hydrophobic (3.1), isomer of leucine (L)', 'P'=>'no charge, aliphatic, promotes turns, cyclic, not hydrophobic (-0.29) ', 'F'=>'no charge, aromatic, absorbs UV, hydrophobic (2.5)', 'Y'=> 'weak charge, aromatic, absorbs UV, hydrogen bonding, not hydrophilic (0.08) ', 'W'=> 'largest amino acid and rarest amino acid, no charge, aromatic, absorbs UV, hydrogen bonding, hydrophobic (1.5)', 'S'=>'no charge, polar, hydrogen bonding, hydrophilic (-1.1)', 'T'=>'no charge, polar, hydrogen bonding, hydrophilic (-0.75)', 'N'=>'amide of Aspartate (D), polar, hydrogen bonding, no charge, hydrophilic (-2.7)', 'Q'=>'amide of Glutamate (E), polar, hydrogen bonding, no charge, hydrophilic (-2.9)', 'C'=>'Sulfur analog of Serine(S), weak charge, forms disulfide bonds, not hydrophilic (0.17)', 'M'=>'initiator of proteins, containing Sulfur, no charge, hydrophobic (1.1)' , 'D'=>'acidic, negative charge, hydrophilic (-3.0)', 'E'=>'acidic, negative charge, hydrophilic (-3.0)', 'H'=>'imidizole in side chain, basic, reactive, weak positive charge, hydrophilic (-1.7)' , 'K'=> 'amine in side chain, basic, reactive, strong positive charge, hydrophilic (-4.6)' , 'R'=> 'guanidinium side chain, basic, strongest positive charge, hydrophilic (-7.5)', '*'=>'stop' , '-'=>'deletion');
    if(exists $p {
        $aa
    }
    ) {
        return $p {
            $aa
        };
    }
    elsif(length($aa)) {
        print STDERR "Bad aa \"$aa\"!!\n";
        exit;
    }
}
sub prot_modif_type {
    my ($sequence_protref,$sequence_prot,$mut_hash,$syn,$chrom,$nhap) = @_;
    my $color="green";
    my $genomique=0;
    my %hash_mut_type;
    my $haplo="no_mutations";
    my $info_haplo="no_mutations";
    my $nb_diff=0;
    my %mut_hash_aa=%$mut_hash;
    for (my $i=0;$i<length($sequence_protref);
    $i++) {
        my $aaref=substr($sequence_protref,$i,1);
        my $aamut=substr($sequence_prot,$i,1);
        my $pos = $i;
        $pos+=1;
        # print $pos."position###\n";
        if (defined $mut_hash_aa {
            $nhap
        } {
            $pos
        }
        ) {
            $genomique=$mut_hash_aa {
                $nhap
            } {
                $pos
            } {
                "vcf"
            };
            # $hash_mut_type{$pos}{"frameshift"}=0;
            $hash_mut_type {
                $pos
            } {
                "genomique"
            }
            =$genomique;
            #  print $mut_hash_aa{$nhap}{$pos}." phase\n";
            if ($mut_hash_aa {
                $nhap
            } {
                $pos
            } {
                "phase"
            }
            eq "missing") {
                $hash_mut_type {
                    $pos
                } {
                    "color"
                }
                ="red";
                $hash_mut_type {
                    $pos
                } {
                    "type"
                }
                ="missing";
                $hash_mut_type {
                    $pos
                } {
                    "frameshift"
                }
                =0;
                $hash_mut_type {
                    $pos
                } {
                    "ref"
                }
                =".";
                $hash_mut_type {
                    $pos
                } {
                    "mut"
                }
                =".";
                $hash_mut_type {
                    $pos
                } {
                    "symbol"
                }
                =5;
                $hash_mut_type {
                    $pos
                } {
                    "genomique"
                }
                =$genomique;
                $nb_diff++;
                next;
            }
            elsif ($nhap eq "hap1" && $mut_hash_aa {
                $nhap
            } {
                $pos
            } {
                "phase"
            }
            eq "yes") {
                $color="blue2";
                $hash_mut_type {
                    $pos
                } {
                    "color"
                }
                ="blue2";
            }
            elsif ($nhap eq "hap2" && $mut_hash_aa {
                $nhap
            } {
                $pos
            } {
                "phase"
            }
            eq "yes") {
                $color="indianred3";
                $hash_mut_type {
                    $pos
                } {
                    "color"
                }
                ="orangered2";
            }
            else {
                $color="azure4";
                $hash_mut_type {
                    $pos
                } {
                    "color"
                }
                ="azure4";
            }
            #print $hash_mut_type{$pos}{"color"}." couleur\n";
        }
        else {
            $hash_mut_type {
                $pos
            } {
                "genomique"
            }
            =$genomique;
            $hash_mut_type {
                $pos
            } {
                "color"
            }
            =$color;
        }
        $hash_mut_type {
            $pos
        } {
            "frameshift"
        }
        =0;
        if ($pos==1 && $aaref ne $aamut) {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            =$aaref;
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            =$aamut;
            $hash_mut_type {
                $pos
            } {
                "type"
            }
            ="start-loss";
            $hash_mut_type {
                $pos
            } {
                "symbol"
            }
            =4;
            $i = length($sequence_protref);
            $nb_diff++;
            next;
        }
        if ($aaref eq $aamut) {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            =$aaref;
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            =$aamut;
            if ($syn eq 'yes' && defined $mut_hash_aa {
                $nhap
            } {
                $pos
            }
            ) {
                #   print $mut_hash_aa{$nhap}{$pos}."\n";
                $hash_mut_type {
                    $pos
                } {
                    "type"
                }
                ="synonymous";
                $hash_mut_type {
                    $pos
                } {
                    "symbol"
                }
                =1;
                $nb_diff++;
            }
            else {
                $hash_mut_type {
                    $pos
                } {
                    "type"
                }
                ="match";
            }
        }
        elsif ($aaref eq '*') {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            ="*";
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            =$aamut;
            $hash_mut_type {
                $pos
            } {
                "type"
            }
            ="elongation";
            $hash_mut_type {
                $pos
            } {
                "symbol"
            }
            =13;
            $i=length($sequence_protref);
            $nb_diff++;
        }
        elsif ($aamut eq '*') {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            =$aaref;
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            ="*";
            $hash_mut_type {
                $pos
            } {
                "type"
            }
            ="stop";
            $hash_mut_type {
                $pos
            } {
                "symbol"
            }
            =8;
            $i=length($sequence_protref);
            $nb_diff++;
        }
        elsif ($aamut eq '-') {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            =$aaref;
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            ="del";
            $hash_mut_type {
                $pos
            } {
                "type"
            }
            ="deletion";
            $hash_mut_type {
                $pos
            } {
                "symbol"
            }
            =17;
            if ($hash_mut_type {
                $pos-1
            } {
                "type"
            }
            eq"synonymous") {
                $hash_mut_type {
                    $pos-1
                } {
                    "type"
                }
                ="match";
            }
            $nb_diff++;
        }
        elsif ($aaref eq '-') {
            if (($pos-1)== 0
                ){
                $hash_mut_type {
                    ($pos-1)
                } {
                    "mut"
                }
                =$hash_mut_type {
                    ($pos-1)
                } {
                    "mut"
                }
                .$aamut;
                $hash_mut_type {
                    ($pos-1)
                } {
                    "type"
                }
                ="insertion";
                $hash_mut_type {
                    $pos-1
                } {
                    "ref"
                }
                ="ins";
                $hash_mut_type {
                    $pos-1
                } {
                    "symbol"
                }
                =6;
                $nb_diff++;
            }
            else {
                $hash_mut_type {
                    ($pos-1)
                } {
                    "mut"
                }
                =$hash_mut_type {
                    ($pos-1)
                } {
                    "mut"
                }
                .$aamut;
                $hash_mut_type {
                    ($pos-1)
                } {
                    "type"
                }
                ="insertion";
                $hash_mut_type {
                    $pos-1
                } {
                    "ref"
                }
                ="ins";
                $hash_mut_type {
                    $pos-1
                } {
                    "symbol"
                }
                =6;
                $nb_diff++;
            }
            substr($sequence_protref,$i,1,"");
            substr($sequence_prot,$i,1,"");
            $i--;
        }
        else {
            $hash_mut_type {
                $pos
            } {
                "ref"
            }
            =$aaref;
            $hash_mut_type {
                $pos
            } {
                "mut"
            }
            =$aamut;
            $hash_mut_type {
                $pos
            } {
                "type"
            }
            ="mismatch";
            $hash_mut_type {
                $pos
            } {
                "symbol"
            }
            =19;
            if (defined $hash_mut_type {
                $pos-1
            } {
                "type"
            }
            ) {
                if ($hash_mut_type {
                    $pos-1
                } {
                    "type"
                }
                eq"mismatch" && !defined $mut_hash_aa {
                    $nhap
                } {
                    $pos
                }
                ) {
                    $hash_mut_type {
                        $pos
                    } {
                        "frameshift"
                    }
                    =$hash_mut_type {
                        $pos-1
                    } {
                        "frameshift"
                    }
                    +1;
                    if ($hash_mut_type {
                        $pos-1
                    } {
                        "frameshift"
                    }
                    ==0) {
                        $hash_mut_type {
                            $pos-1
                        } {
                            "frameshift"
                        }
                        =-1;
                    }
                }
                elsif ($hash_mut_type {
                    $pos-2
                } {
                    "type"
                }
                eq"mismatch" && $hash_mut_type {
                    $pos-2
                } {
                    "type"
                }
                ne"synonymous" && !defined $mut_hash_aa {
                    $nhap
                } {
                    $pos
                }
                ) {
                    $hash_mut_type {
                        $pos-1
                    } {
                        "frameshift"
                    }
                    =$hash_mut_type {
                        $pos-2
                    } {
                        "frameshift"
                    }
                    +1;
                    $hash_mut_type {
                        $pos
                    } {
                        "frameshift"
                    }
                    =$hash_mut_type {
                        $pos-2
                    } {
                        "frameshift"
                    }
                    +2;
                    if ($hash_mut_type {
                        $pos-2
                    } {
                        "frameshift"
                    }
                    ==0) {
                        $hash_mut_type {
                            $pos-2
                        } {
                            "frameshift"
                        }
                        =-1;
                    }
                }
                $nb_diff++;
            }
            else {
                delete $hash_mut_type {
                    $pos-1
                };
            }
        }
    }
    #print Dumper (\%hash_mut_type);
	#print Dumper ($hash_mut_type{1});
	#print Dumper ($hash_mut_type{0});
    if ($nb_diff!= 0
        ){
        $haplo="";
        $info_haplo="";
        my $prev_k;
        my $fsize=0;
        for my $k( sort {
            $a<=>$b
        }
        keys %hash_mut_type) {
            if ($hash_mut_type {
                $k
            } {
                "type"
            }
            ne "match") {
                #print $k." position mutation \n".$hash_mut_type{$k}{"type"}." type \n".$hash_mut_type{$k}{"frameshift"}." frameshift\n";
                if ($hash_mut_type {
                    $k
                } {
                    "frameshift"
                }
                !=0 && $hash_mut_type {
                    $k
                } {
                    "frameshift"
                }
                !=-1) {
                    $fsize=$hash_mut_type {
                        $k
                    } {
                        "frameshift"
                    };
                }
                if ($hash_mut_type {
                    $k
                } {
                    "frameshift"
                }
                ==-1 || $hash_mut_type {
                    $k
                } {
                    "frameshift"
                }
                ==0) {
                    if ($fsize!=0) {
                        if ($haplo eq "" && $prev_k != 0
                            ){
                            $info_haplo .= "FS".$fsize.":".($prev_k-$fsize).":".$prev_k.":frameshift:19:".$hash_mut_type {
                                $k
                            } {
                                "color"
                            };
                            $haplo .= "FS".$fsize.":".($prev_k-$fsize).":".$prev_k;
                            $fsize=0;
                        }
                        elsif ($prev_k !=0) {
                            $info_haplo .= ",FS".$fsize.":".($prev_k-$fsize).":".$prev_k.":frameshift:19:".$hash_mut_type {
                                $k
                            } {
                                "color"
                            };
                            $haplo .= ",FS".$fsize.":".($prev_k-$fsize).":".$prev_k;
                            $fsize=0;
                        }
                    }
                    if ($haplo eq "" && $hash_mut_type {
                        $k
                    } {
                        "frameshift"
                    }
                    !=-1) {
                        $info_haplo = $hash_mut_type {
                            $k
                        } {
                            "ref"
                        }
                        .":".$k.":".$hash_mut_type {
                            $k
                        } {
                            "mut"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "type"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "symbol"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "color"
                        };
                        $haplo = $hash_mut_type {
                            $k
                        } {
                            "ref"
                        }
                        .":".$k.":".$hash_mut_type {
                            $k
                        } {
                            "mut"
                        };
                    }
                    elsif ($hash_mut_type {
                        $k
                    } {
                        "frameshift"
                    }
                    !=-1) {
                        $info_haplo .= ",".$hash_mut_type {
                            $k
                        } {
                            "ref"
                        }
                        .":".$k.":".$hash_mut_type {
                            $k
                        } {
                            "mut"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "type"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "symbol"
                        }
                        .":".$hash_mut_type {
                            $k
                        } {
                            "color"
                        };
                        $haplo .= ",".$hash_mut_type {
                            $k
                        } {
                            "ref"
                        }
                        .":".$k.":".$hash_mut_type {
                            $k
                        } {
                            "mut"
                        };
                    }
                }
                $prev_k=$k;
            }
        }
    }
    #tludwig : do the following code do anything ?
    for my $m (sort {
        $a<=>$b
    }
    keys (%hash_mut_type)) {
        if ($hash_mut_type {
            $m
        } {
            "type"
        }
        ne "match") {
            #print $m." ".Dumper(\$hash_mut_type{$m})."\n";
        }
    }
    return ($info_haplo, $haplo, \%hash_mut_type);
}
sub html_file {
    my ($chr, $aa_pos, $htlm_file, $pdf, $nb1, $nb2, $syn, $mut_hash_ref, $hash_mut_type1_ref, $hash_mut_type2_ref, $clin_hash_ref) = @_;
    my %mut_hash=%$mut_hash_ref;
    my %hash_mut_type1=%$hash_mut_type1_ref;
    my %hash_mut_type2=%$hash_mut_type2_ref;
    my %clin_hash=%$clin_hash_ref;
    #print Dumper(\%mut_hash);
    open HTML, ">>".$htlm_file;
    print HTML "<p class=\"haplomodification\"><h3>Haplotype 1 modifications</h3> \n";
    if ($nb1 ==0) {
        print HTML "<br><blockquote>No changes from reference sequence\n</blockquote><br>";
    }
    else {
        print HTML "<ul class=\"ulmut\">";
        my $prev_gen=0;
        my $stop_gen_pos="no";
        for my $m (sort {
            $a<=>$b
        }
        keys % {
            hash_mut_type1
        }
        ) {
            if ($hash_mut_type1 {
                $m
            } {
                "type"
            }
            ne "match") {
                if ($hash_mut_type1 {
                    $m
                } {
                    "genomique"
                }
                !=$prev_gen) {
                    print HTML "<br><li>".$mut_hash {
                        $chr
                    } {
                        $hash_mut_type1 {
                            $m
                        } {
                            "genomique"
                        }
                    } {
                        'hap1'
                    } {
                        'vcf'
                    }
                    ." amino acid position ".$m." \n</li><br><blockquote>";
                    $stop_gen_pos="no";
                }
                my $frame= $mut_hash {
                    $chr
                } {
                    $hash_mut_type1 {
                        $m
                    } {
                        "genomique"
                    }
                } {
                    'hap1'
                } {
                    'amino'
                };
                if($stop_gen_pos eq "yes") {
                    next;
                }
                #	print HTML "<br><li>hap 1 ".$chr.":".$chr.$m.":".$mut_hash{$chr}{$m}{'hap1'}{'nt'}." amino acid position ".$mut_hash{$chr}{$m}{'hap1'}{'amino'}." \n</li><br><blockquote>";
                elsif ($hash_mut_type1 {
                    $m
                } {
                    "frameshift"
                }
                == -1) {
                    $stop_gen_pos="yes";
                    while ($hash_mut_type1 {
                        $frame
                    } {
                        "type"
                    }
                    ne "elongation" && $hash_mut_type1 {
                        $frame
                    } {
                        "type"
                    }
                    ne "stop" && $frame != $aa_pos) {
                        $frame++;
                    }
                    if ($hash_mut_type1 {
                        $frame
                    } {
                        "type"
                    }
                    eq "elongation") {
                        #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                        print HTML "Frameshift mutation that results in protein elongation <br>";
                    }
                    elsif ($hash_mut_type1 {
                        $frame
                    } {
                        "type"
                    }
                    eq "stop") {
                        #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                        print HTML "Frameshift mutation that results in a premature stop codon at position ".$frame."<br>";
                    }
                    elsif ($frame == $aa_pos) {
                        print HTML "Frameshift mutation with no elongation of protein or premature stop codon <br>";
                    }
                }
                elsif ($hash_mut_type1 {
                    $m
                } {
                    "type"
                }
                eq "missing") {
                    print HTML "Missing Data<br>";
                }
                elsif ($hash_mut_type1 {
                    $m
                } {
                    "type"
                }
                eq "mismatch") {
                    my $aapref=aa2properties($hash_mut_type1 {
                        $m
                    } {
                        "ref"
                    }
                    );
                    my $aapmut=aa2properties($hash_mut_type1 {
                        $m
                    } {
                        "mut"
                    }
                    );
                    #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                    print HTML $hash_mut_type1 {
                        $m
                    } {
                        "ref"
                    }
                    ."=>".$hash_mut_type1 {
                        $m
                    } {
                        "mut"
                    }
                    ." ".$hash_mut_type1 {
                        $m
                    } {
                        "type"
                    }
                    ." mutation<br>";
                    print HTML "Change from ".$aapref."amino acid to ".$aapmut."<br>";
                }
                else {
                    # print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                    print HTML $hash_mut_type1 {
                        $m
                    } {
                        "ref"
                    }
                    ."=>".$hash_mut_type1 {
                        $m
                    } {
                        "mut"
                    }
                    ." ".$hash_mut_type1 {
                        $m
                    } {
                        "type"
                    }
                    ." mutation<br>";
                }
                if(exists( $clin_hash {
                    $chr
                } {
                    $hash_mut_type1 {
                        $m
                    } {
                        "genomique"
                    }
                }
                )) {
                    print HTML " ClinVar Information at this position : \n<br>";
                    print HTML "<blockquote>";
                    foreach my $clin_info (keys % {
                        $clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        }
                    }
                    ) {
                        print HTML "<ul class=\"ulclinvar\">";
                        #if ($clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{$clin_info}{"ALT"} =~ /HGVS/){
                        print HTML "<li>".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "RS"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "NC"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "CLNSIG"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "CLNDBN"
                        }
                        ." <a href=\"https://www.ncbi.nlm.nih.gov/clinvar/variation/".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type1 {
                                $m
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "ID"
                        }
                        ."\">"."see clinvar"."</a> </li><br>";
                        #print HTML "<li>".$clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{$clin_info}{"RS"}." <a href=\"https://www.ncbi.nlm.nih.gov/clinvar/variation/".$clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{"ID"}."\">".$clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{$clin_info}{"NC"}."</a>  ".$clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{$clin_info}{"CLNSIG"}."   ".$clin_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{$clin_info}{"CLNDBN"}."</li><br>";
						#}
                        print HTML "</ul>";
                    }
                    print HTML "</blockquote>";
                }
                print HTML "</blockquote>";
                $prev_gen=$hash_mut_type1 {
                    $m
                } {
                    "genomique"
                };
            }
        }
        print HTML "</ul>"
    }
    print HTML "</p>";
    print HTML "<p class=\"haplomodification\"><h3>Haplotype 2 modifications </h3> \n";
    if ($nb2 ==0) {
        print HTML "<br><blockquote>No changes from reference sequence\n</blockquote><br>";
    }
    else {
        print HTML "<ul class=\"ulmut\">";
        my $prev_gen2=0;
        my $stop_gen_pos2="no";
        for my $m2 (sort {
            $a<=>$b
        }
        keys % {
            hash_mut_type2
        }
        ) {
            if ($hash_mut_type2 {
                $m2
            } {
                "type"
            }
            ne "match") {
                if ($hash_mut_type2 {
                    $m2
                } {
                    "genomique"
                }
                !=$prev_gen2) {
                    print HTML "<br><li>".$mut_hash {
                        $chr
                    } {
                        $hash_mut_type2 {
                            $m2
                        } {
                            "genomique"
                        }
                    } {
                        'hap2'
                    } {
                        'vcf'
                    }
                    ." amino acid position ".$m2." \n</li><br><blockquote>";
                    $stop_gen_pos2="no";
                }
                my $frame2= $mut_hash {
                    $chr
                } {
                    $hash_mut_type2 {
                        $m2
                    } {
                        "genomique"
                    }
                } {
                    'hap2'
                } {
                    'amino'
                };
                if($stop_gen_pos2 eq "yes") {
                    next;
                }
                #print HTML "<br><li>hap 2 ".$chr.":".$chr.$m2.":".$mut_hash{$chr}{$m2}{'hap2'}{'nt'}." amino acid position ".$mut_hash{$chr}{$m2}{'hap2'}{'amino'}." \n</li><br><blockquote>";
				#print $hash_mut_type2{$mut_hash{$chr}{$m2}{'hap2'}{'amino'}}{"ref"}.">".$hash_mut_type2{$mut_hash{$chr}{$m2}{'hap2'}{'amino'}}{"mut"}." ".$hash_mut_type2{$mut_hash{$chr}{$m2}{#'hap2'}{'amino'}}{"type"}." mutation<br>\n";
                elsif ($hash_mut_type2 {
                    $m2
                } {
                    "frameshift"
                }
                == -1) {
                    $stop_gen_pos2="yes";
                    while ($hash_mut_type2 {
                        $frame2
                    } {
                        "type"
                    }
                    ne "elongation" && $hash_mut_type2 {
                        $frame2
                    } {
                        "type"
                    }
                    ne "stop" && $frame2 != $aa_pos) {
                        $frame2++;
                    }
                    if ($hash_mut_type2 {
                        $frame2
                    } {
                        "type"
                    }
                    eq "elongation") {
                        #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                        print HTML "Frameshift mutation that results in protein elongation <br>";
                    }
                    elsif ($hash_mut_type2 {
                        $frame2
                    } {
                        "type"
                    }
                    eq "stop") {
                        #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                        print HTML "Frameshift mutation that results in a premature stop codon at position ".$frame2."<br>";
                    }
                    elsif ($frame2 == $aa_pos) {
                        print HTML "Frameshift mutation with no elongation of protein or premature stop codon <br>";
                    }
                }
                elsif ($hash_mut_type2 {
                    $m2
                } {
                    "type"
                }
                eq "missing") {
                    print HTML "Missing Data<br>";
                }
                elsif ($hash_mut_type2 {
                    $m2
                } {
                    "type"
                }
                eq "mismatch") {
                    my $aapref=aa2properties($hash_mut_type2 {
                        $m2
                    } {
                        "ref"
                    }
                    );
                    my $aapmut=aa2properties($hash_mut_type2 {
                        $m2
                    } {
                        "mut"
                    }
                    );
                    #print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                    print HTML $hash_mut_type2 {
                        $m2
                    } {
                        "ref"
                    }
                    ."=>".$hash_mut_type2 {
                        $m2
                    } {
                        "mut"
                    }
                    ." ".$hash_mut_type2 {
                        $m2
                    } {
                        "type"
                    }
                    ." mutation<br>";
                    print HTML "Change from ".$aapref."amino acid to ".$aapmut."<br>";
                }
                else {
                    # print HTML "<br><li>".$mut_hash{$chr}{$hash_mut_type1{$m}{"genomique"}}{'hap1'}{'vcf'}." amino acid position ".$m." \n</li><br><blockquote>";
                    print HTML $hash_mut_type2 {
                        $m2
                    } {
                        "ref"
                    }
                    ."=>".$hash_mut_type2 {
                        $m2
                    } {
                        "mut"
                    }
                    ." ".$hash_mut_type2 {
                        $m2
                    } {
                        "type"
                    }
                    ." mutation<br>";
                }
                if(exists( $clin_hash {
                    $chr
                } {
                    $hash_mut_type2 {
                        $m2
                    } {
                        "genomique"
                    }
                }
                )) {
                    print HTML " ClinVar Information at this position : \n<br>";
                    print HTML "<blockquote>";
                    foreach my $clin_info (keys % {
                        $clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        }
                    }
                    ) {
                        print HTML "<ul class=\"ulclinvar\">";
                        #if ($clin_info =~ /HGVS/){
                        print HTML "<li>".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "RS"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "NC"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "CLNSIG"
                        }
                        ."   ".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "CLNDBN"
                        }
                        ." <a href=\"https://www.ncbi.nlm.nih.gov/clinvar/variation/".$clin_hash {
                            $chr
                        } {
                            $hash_mut_type2 {
                                $m2
                            } {
                                "genomique"
                            }
                        } {
                            $clin_info
                        } {
                            "ID"
                        }
                        ."\">"."see clinvar"."</a>  </li><br>";
                        #print HTML "<li>".$clin_hash{$chr}{$hash_mut_type2{$m2}{"genomique"}}{$clin_info}{"RS"}." <a href=\"https://www.ncbi.nlm.nih.gov/clinvar/variation/".$clin_hash{$chr}{$hash_mut_type2{$m2}{"genomique"}}{"ID"}."\">".$clin_hash{$chr}{$hash_mut_type2{$m2}{"genomique"}}{$clin_info}{"NC"}."</a>  ".$clin_hash{$chr}{$hash_mut_type2{$m2}{"genomique"}}{$clin_info}{"CLNSIG"}."   ".$clin_hash{$chr}{$hash_mut_type2{$m2}{"genomique"}}{$clin_info}{"CLNDBN"}."</li><br>";
						#}
                        print HTML "</ul>";
                    }
                    print HTML "</blockquote>";
                }
                print HTML "</blockquote>";
                $prev_gen2=$hash_mut_type2 {
                    $m2
                } {
                    "genomique"
                };
            }
        }
        print HTML "</ul>"
    }
    print HTML "</p>";
    print HTML "<br><object class=\"bigpdf\" data=\"".$pdf."\" type=\"application/pdf\">alt : <a href=\"".$pdf."\">protein.pdf</a></object>";
    print HTML "<div w3-include-html=\"http://lysine.univ-brest.fr/GEMPROT/footer.html\" id=\"gemprotfooter\" class=\"gemprotfooter\"></div>";
    print HTML "</p>\n<script>includeHTML();</script></body>\n</html>";
    close HTML;
}
# Get pfam result to hash table
############
sub get_url_pfam {
    my ($sequence_gene_prot, $gene, $cache) = @_;
    my @url;
    my @result;
    my @empty;
    my %prot_domain;
    my $start =0;
    my $pfam_store_file = $cache."/pfam_store.txt";
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
    my $date_release = $root->getAttribute( 'release_date' );
    chomp $date_release;
    if (!-e $pfam_store_file) {
        #no pfam_store file
        open PSF, ">>".$pfam_store_file;
        #tludwig open file only here
		#printf "Existe pas Release_date= ".$date_release ."\n";
        print PSF "Release_date=".$date_release."\n";
        print PSF "Gene\tDomaine\tStart\tEnd\tAccession\n";
        close PSF;
        #tludwig close file
    }
    else {
        open(FL, "grep \"Release\" ".$pfam_store_file." |");
        my @firstline =<FL>;
        close FL;
        #tludwig close file after reading
		#printf "premiere ligne ".Dumper (\@firstline)."\n";
        my @fline=split(/=/, $firstline[0]);
        #printf Dumper (\@fline);
        my $date_pfam=$fline[1];
        chomp $date_pfam;
        if ($date_pfam ne $date_release) {
            #printf 'Existe mais pas bon Release_date= '.$date_release ."\n";
            system("rm ".$pfam_store_file);
            #printf "Release_date= ".$date_release ."\n";
            open PSF, ">>".$pfam_store_file;
            print PSF "Release_date=".$date_release."\n";
            print PSF "Gene\tDomaine\tStart\tEnd\tAccession\n";
            close PSF;
            #tludwig close file
        }
    }
    printf "start pfam\n";
    open(GGF, "grep -w \"".$gene."\" ".$pfam_store_file." |");
    my @genedomain =<GGF>;
    if (@genedomain== 0
        ){
        while (@result == 0
            ){
            @url=`curl --silent -LH 'Expect:' -F seq='<$sequence_gene_prot' -F output=xml 'http://pfam.xfam.org/search/sequence'`
            ;
			#printf "curl --silent -LH 'Expect:' -F seq='<$sequence_gene_prot' -F output=xml 'http://pfam.xfam.org/search/sequence'\n";
            @result = grep(/result_url/, @url);
            @empty = grep(/empty/, @url);
            if (@empty != 0
                ){
                $prot_domain {
                    $start
                } {
                    "Start"
                }
                = 1;
                $prot_domain {
                    $start
                } {
                    "End"
                }
                = 1;
                $prot_domain {
                    $start
                } {
                    "Domain"
                }
                = "No_domains";
                $prot_domain {
                    $start
                } {
                    "Accession"
                }
                = 1;
                #print Dumper (\%prot_domain);
                open PSF, ">>".$pfam_store_file;
                print PSF $gene."\t"."No_domains"."\t"."1"."\t"."1"."\t"."1"."\n";
                close PSF;
                return \%prot_domain;
            }
            #printf Dumper (\@url);
			#printf Dumper (\@result);
        }
        $result[0] =~ s/result_url|\>|\<| //g;
        $result[0] =~ s/\/$//;
        my $res_url= $result[0
        ];
        chomp $res_url;
        #$res_url = "'".$res_url."'";
        my $pfam="";
        printf "\n==> RUN PFAM\n";
        while ($pfam !~ m/^\</) {
            $pfam=`curl --silent -LH 'Expect:' -F output=xml 'http://pfam.xfam.org$res_url'`
            ;
			# printf "curl --silent -LH 'Expect:' -F output=xml 'http://pfam.xfam.org$res_url' ";
            if ($pfam !~ m/^\</) {
                print $pfam."   \r";
            }
        }
        print "\n==> END PFAM\n";
        #print "##".$pfam."###";
        if ($pfam !~ /matches/) {
            $prot_domain {
                $start
            } {
                "Start"
            }
            = 1;
            $prot_domain {
                $start
            } {
                "End"
            }
            = 1;
            $prot_domain {
                $start
            } {
                "Domain"
            }
            = "No_domains";
            $prot_domain {
                $start
            } {
                "Accession"
            }
            = 1;
            #print Dumper (\%prot_domain);
            open PSF, ">>".$pfam_store_file;
            print PSF $gene."\t"."No_domains"."\t"."1"."\t"."1"."\t"."1"."\n";
            close PSF;
            return \%prot_domain;
        }
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
            my $ID = 0;
			my $start =0;
			#print Dumper \@entry4;
            my $domain_gene="";
            foreach my $match (@entry4) {
                #print 'domain: ' . $match->getAttribute( 'id' ) . "\n";
                my ( @entry5 ) = $match->getChildrenByTagName( 'location' );
                foreach my $loc (@entry5) {
                    $start=$loc->getAttribute( 'ali_start' );
                    $prot_domain {
                        $start
                    } {
                        "Start"
                    }
                    = $loc->getAttribute( 'ali_start' );
                    $prot_domain {
                        $start
                    } {
                        "End"
                    }
                    = $loc->getAttribute( 'ali_end' );
                    $prot_domain {
                        $start
                    } {
                        "Domain"
                    }
                    = $match->getAttribute( 'id' );
                    $prot_domain {
                        $start
                    } {
                        "Accession"
                    }
                    = $match->getAttribute( 'accession' );
                    my $start=$loc->getAttribute( 'ali_start' );
                    my $end = $loc->getAttribute( 'ali_end' );
                    my $domain = $match->getAttribute( 'id' );
                    my $accession = $match->getAttribute( 'accession' );
                    $ID++;
                    $domain_gene=$domain_gene.$gene."\t".$domain."\t".$start."\t".$end."\t".$accession."\n";
                }
            }
            open PSF, ">>".$pfam_store_file;
            print PSF $domain_gene;
        }
    }
    else {
        foreach my $dom (@genedomain) {
            my @d_line= split (/\t/ , $dom);
            my $start = $d_line[2];
            $prot_domain {
                $start
            } {
                "Start"
            }
            = $d_line[2];
            $prot_domain {
                $start
            } {
                "End"
            }
            = $d_line[3];
            $prot_domain {
                $start
            } {
                "Domain"
            }
            = $d_line[1];
            $prot_domain {
                $start
            } {
                "Accession"
            }
            = $d_line[4];
        }
    }
    close PSF;
    #print Dumper \%prot_domain;
    return \%prot_domain;
}
sub domain_file {
    my ($d_file,$domain,$gene)=@_;
    #printf $domain;
    my %prot_domain =%$domain;
    #print Dumper(\%prot_domain);
    open DFILE, ">".$d_file;
    print DFILE "gene\tdomain_name\tstart_site\tend_site\n";
    #print "domain_name\tstart_site\tend_site\n";
    foreach my $key (sort {
        $a <=> $b
    }
    keys %prot_domain) {
        print DFILE $gene."\t".$prot_domain {
            $key
        } {
            "Domain"
        }
        ."\t".$prot_domain {
            $key
        } {
            "Start"
        }
        ."\t".$prot_domain {
            $key
        } {
            "End"
        }
        ."\n";
    }
    close DFILE;
}
# Get column in vcf header
############
sub get_header {
    my ($line1,$sep,$col_name) = @_;
    chomp($line1);
    my @elmnts=split(/$sep/,$line1);
    my %header;
    my $i=0;
    while ($i<scalar(@elmnts)) {
        $header {
            $elmnts[$i]
        }
        =$i;
        $i++;
    }
    if (defined $header {
        $col_name
    }
    ) {
        my $num_col=$header {
            $col_name
        };
        return $num_col;
    }
    else {
        my $num_col=0;
        return $num_col;
    }
}
# Get chromosome format (chr1 or 1)
###########
sub get_chr {
    my ($vcf) = @_;
    my $with_chr="no";
    if ($vcf =~ /.gz$/) {
        open(IN, "zcat $vcf | head -2000 | grep -v \"^#\" | head -1 | cut -f1 | ") || die "can not open pipe to $vcf";
    }
    else {
        open(IN, "cat $vcf | head -2000 | grep -v \"^#\" | head -1 | cut -f1 | ") || die "can not open pipe to $vcf";
    }
    my $header = <IN>;
    if ($header=~"^chr") {
        $with_chr="yes";
    }
    return $with_chr;
}
# Get sample name when sample-file is not provied
###########
sub get_sample {
    my ($vcf, $sample_rec) = @_;
    my @elmnts;
    if ($vcf =~ /.gz$/) {
        open(IN, "zcat $vcf | head -2000 | grep \"^#\" | tail -1 | ") || die "can not open pipe to $vcf";
    }
    else {
        open(IN, "cat $vcf | head -2000 | grep \"^#\" | tail -1 | ") || die "can not open pipe to $vcf";
    }
    my $header = <IN>;
    @elmnts=split(/\t/,$header);
    my $i=9;
    open SFILE, ">".$sample_rec;
    while ($i<scalar(@elmnts)) {
        if ($i<scalar(@elmnts)-1) {
            print SFILE $elmnts[$i]."\n";
            #printf $elmnts[$i]."\n";
        }
        else {
            print SFILE $elmnts[$i];
        }
        $i++;
    }
}
# Change sequence with variation in vcf
############
sub parsevcf2 {
    my ($vcf,$seq_hash,$seq_hash1,$seq_hash2,$seq_ref1,$seq_ref2,$sample,$clin,$gene) = @_;
    if ($vcf =~ /.gz$/) {
        open(IN, "gunzip -c $vcf |") || die "can’t open pipe to $vcf";
    }
    else {
        open(IN, $vcf) || die "can’t open $vcf";
    }
    #open VCF, $vcf;
    my @variant = <IN>;
    my %seq_hash_ref=%$seq_hash;
    my %seq_hash_hap1=%$seq_hash1;
    my %seq_hash_hap2=%$seq_hash2;
    my %seq_hash_ref1=%$seq_ref1;
    my %seq_hash_ref2=%$seq_ref2;
    my %clin_hash=%$clin;
    my %mut_hash=();
    my %mut_hash_pos=();
    #my @variant = <VCF> ;
	my $mut_nb1=0;
	my $mut_nb2=0;
    my $num_col="";
    my $phased="";
    #print "##########Summary of nucleotide mutations present in $gene sequence\n";
    foreach my $line (@variant) {
        ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  B00FWWD
		#print $line." ligne du variant \n\n";
		next if ( $line =~ /(^\s*$)|(^\")|(^-)|(^:)|(^#)/ && $line!~"CHROM");
        chomp($line);
        my @fich1 = split(/\t/, $line);
        my $variant_pos=$fich1[1];
        if ($line=~"^#CHROM") {
            $num_col=get_header($line,'\\t',$sample);
            if ($num_col== 0
                ){
                printf "##########ERROR MESSAGE : sample $sample doesn't exist in the vcf\n";
                return (\%seq_hash_hap1, \%seq_hash_hap2, \%seq_hash_ref1, \%seq_hash_ref2, \%mut_hash, \%mut_hash_pos, $mut_nb1, $mut_nb2);
            }
        }
        elsif(defined $seq_hash_ref {
            $fich1[
        0
        ]} {
            $variant_pos
        }
        ) {
            # print $fich1[3]." ".$fich1[4]." position $fich1[1]in coding region\n ";
			# regex for ALT and REF contain only nucleotide (ATGC) and add "," for billelic multiallelic site
            if ($fich1[3]=~ m/^[ATGC\,]+$/ && $fich1[4]=~ m/^[ATGC\,]+$/) {
                my @geno = split(/:/, $fich1[$num_col]);
                # verify if the genotype is phase or not
                if ($geno[0] =~ m/\|/) {
                    $phased="yes";
                }
                else {
                    $phased="no";
                    $fich1[$num_col]=~ tr/\//\|/;
                }
                #print "ce fichier est phase : $phased \n";
				#print $fich1[$num_col];
				# print $fich1[3]." ".$fich1[4]." position $fich1[1] is ATGC variation\n";
                my @info = split(/[:|]/, $fich1[$num_col]);
                my $i=0;
				my $i_inf=0; # last nt before ALT > REF => start insertion
                my $geno;
                if($geno[0] =~ m/\./) {
                    $mut_nb1++;
                    $mut_nb2++;
                    $mut_hash_pos {
                        "hap1"
                    } {
                        $seq_hash_ref {
                            $fich1[
                        0
                        ]} {
                            $variant_pos+$i
                        } {
                            "amino"
                        }
                    } {
                        "phase"
                    }
                    = "missing";
                    $mut_hash_pos {
                        "hap1"
                    } {
                        $seq_hash_ref {
                            $fich1[
                        0
                        ]} {
                            $variant_pos+$i
                        } {
                            "amino"
                        }
                    } {
                        "vcf"
                    }
                    = $fich1[1];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap1"
                    } {
                        "vcf"
                    }
                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$fich1[4];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap1"
                    } {
                        "nt"
                    }
                    = $fich1[3].">".$fich1[4];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap1"
                    } {
                        "amino"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "amino"
                    };
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap1"
                    } {
                        "exon"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "exon"
                    };
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap1"
                    } {
                        "nt_pos"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "nt_pos"
                    };
                    $mut_hash_pos {
                        "hap2"
                    } {
                        $seq_hash_ref {
                            $fich1[
                        0
                        ]} {
                            $variant_pos+$i
                        } {
                            "amino"
                        }
                    } {
                        "phase"
                    }
                    = "missing";
                    $mut_hash_pos {
                        "hap2"
                    } {
                        $seq_hash_ref {
                            $fich1[
                        0
                        ]} {
                            $variant_pos+$i
                        } {
                            "amino"
                        }
                    } {
                        "vcf"
                    }
                    = $fich1[1];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap2"
                    } {
                        "vcf"
                    }
                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$fich1[4];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap2"
                    } {
                        "nt"
                    }
                    = $fich1[3].">".$fich1[4];
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap2"
                    } {
                        "amino"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "amino"
                    };
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap2"
                    } {
                        "exon"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "exon"
                    };
                    $mut_hash {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "hap2"
                    } {
                        "nt_pos"
                    }
                    = $seq_hash_ref {
                        $fich1[
                    0
                    ]} {
                        $variant_pos+$i
                    } {
                        "nt_pos"
                    };
                }
                else {
                    # Keep genotype for biallelic and multiallelic site
                    if($info[0]!= 0
                        ){
                        $mut_nb1++;
                        #print Dumper (\@info);
						#print "genotype hap1 \n".$info[0]."\n";
                        my @ALT_list = split(/,/, $fich1[4]);
                        #print "\nALT ".$fich1[4]."\n";
						#print Dumper (\@ALT_list);
                        $geno=$info[0]-1;
                        my $ALT = $ALT_list[$geno];
                        #print "\nALT ".$ALT."ALT\n";
                        my @ALT_nt = split(//, $ALT);
                        #print Dumper (\@ALT_nt);
                        my @REF_nt = split(//, $fich1[3]);
                        #print Dumper (\@REF_nt);
						#print length($fich1[3])." taille REF ".$fich1[3]."\n";
                        if (length($ALT)>= length($fich1[3])) {
                            #print $i." i ".(length($fich1[3])-1)."\n";
                            for ($i=0; $i<length($ALT);
                            $i++) {
                                #print $i." i ".(length($fich1[3])-1)."\n";
                                if ($i<=((length($fich1[3]))-1)) {
                                    my $i_inf=$i;
                                    #print "\n c'est un SNP ";
									#if ($REF_nt[$i] ne $ALT_nt[$i]){
                                    $seq_hash_hap1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    }
                                    = $ALT_nt[$i];
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[1];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].">".$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                    #}
                                }
                                else {
                                    if(!defined $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    ) {
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap1"
                                        } {
                                            "vcf"
                                        }
                                        = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap1"
                                        } {
                                            "nt"
                                        }
                                        = $REF_nt[$i_inf].">".$ALT_nt[$i_inf];
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap1"
                                        } {
                                            "amino"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "amino"
                                        };
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap1"
                                        } {
                                            "exon"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "exon"
                                        };
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap1"
                                        } {
                                            "nt_pos"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "nt_pos"
                                        };
                                    }
                                    #print "\n c'est une insertion ";
                                    $seq_hash_ref1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    .="-";
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[1];
                                    $seq_hash_hap1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    = $seq_hash_hap1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    .$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    = $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    .$ALT_nt[$i];
                                }
                            }
                        }
                        else {
                            #print "\n c'est une délétion ";
                            for ($i=0; $i<length($fich1[3]);
                            $i++) {
                                if ($i<=((length($ALT))-1)) {
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].">".$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                }
                                #if ($i>=((length($ALT))-1)) {
								#if (!exists $ALT_nt[$i]){
                                else {
                                    $seq_hash_hap1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    }
                                    ="-";
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $mut_hash_pos {
                                        "hap1"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[1];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].">-";
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap1"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                }
                                #}
								#}
                            }
                        }
                    }
                    if($info[1]!= 0
                        ){
                        $mut_nb2++;
                        #print Dumper (\@info);
						#print "genotype hap2 \n".$info[1]."\n";
                        my @ALT_list = split(/,/, $fich1[4]);
                        #print "\nALT ".$fich1[4]."\n";
						#print Dumper (\@ALT_list);
                        $geno=$info[1]-1;
                        my $ALT = $ALT_list[$geno];
                        #print "\nALT ".$ALT."ALT\n";
                        my @ALT_nt = split(//, $ALT);
                        #	print Dumper (\@ALT_nt);
                        my @REF_nt = split(//, $fich1[3]);
                        #	print Dumper (\@REF_nt);
						#print length($fich1[3])." taille REF ".$fich1[3]."\n";
                        if (length($ALT)>= length($fich1[3])) {
                            #print $i." i ".(length($fich1[3])-1)."\n";
                            for ($i=0; $i<length($ALT);
                            $i++) {
                                #print $i." i ".(length($fich1[3])-1)."\n";
                                if ($i<=(length($fich1[3])-1)) {
                                    my $i_inf=$i;
                                    #	print "\n c'est un SNP ";
									#if ($REF_nt[$i] ne $ALT_nt[$i]){
                                    $seq_hash_hap2 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    }
                                    =$ALT_nt[$i];
                                    $mut_hash_pos {
                                        "hap2"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $mut_hash_pos {
                                        "hap2"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[1];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].">".$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                    #}
                                }
                                else {
                                    if(!defined $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    ) {
                                        $mut_hash_pos {
                                            "hap2"
                                        } {
                                            $seq_hash_ref {
                                                $fich1[
                                            0
                                            ]} {
                                                $variant_pos+$i_inf
                                            } {
                                                "amino"
                                            }
                                        } {
                                            "phase"
                                        }
                                        = $phased;
                                        $mut_hash_pos {
                                            "hap2"
                                        } {
                                            $seq_hash_ref {
                                                $fich1[
                                            0
                                            ]} {
                                                $variant_pos+$i_inf
                                            } {
                                                "amino"
                                            }
                                        } {
                                            "vcf"
                                        }
                                        = $fich1[1];
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap2"
                                        } {
                                            "vcf"
                                        }
                                        = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap2"
                                        } {
                                            "nt"
                                        }
                                        = $REF_nt[$i_inf].">".$ALT_nt[$i_inf];
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap2"
                                        } {
                                            "amino"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "amino"
                                        };
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap2"
                                        } {
                                            "exon"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "exon"
                                        };
                                        $mut_hash {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "hap2"
                                        } {
                                            "nt_pos"
                                        }
                                        = $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "nt_pos"
                                        };
                                    }
                                    #print "\n c'est une insertion ";
                                    $seq_hash_ref2 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    .="-";
                                    $mut_hash_pos {
                                        "hap2"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i_inf
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $seq_hash_hap2 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    = $seq_hash_hap1 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    }
                                    .$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    = $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i_inf
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    .$ALT_nt[$i];
                                }
                            }
                        }
                        else {
                            #print "\n c'est une délétion ";
                            for ($i=0; $i<length($fich1[3]);
                            $i++) {
                                if ($i<=((length($ALT))-1)) {
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].$REF_nt[$i].">".$ALT_nt[$i];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                }
                                #if ($i>=((length($ALT))-1)) {
								#	if (!exists $ALT_nt[$i]){
                                else {
                                    $seq_hash_hap2 {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    }
                                    ="-";
                                    $mut_hash_pos {
                                        "hap2"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "phase"
                                    }
                                    = $phased;
                                    $mut_hash_pos {
                                        "hap2"
                                    } {
                                        $seq_hash_ref {
                                            $fich1[
                                        0
                                        ]} {
                                            $variant_pos+$i
                                        } {
                                            "amino"
                                        }
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[1];
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "vcf"
                                    }
                                    = $fich1[0].":".$fich1[1]." ".$fich1[3].">".$ALT;
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt"
                                    }
                                    = $REF_nt[$i].">-";
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "amino"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "amino"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "exon"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "exon"
                                    };
                                    $mut_hash {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "hap2"
                                    } {
                                        "nt_pos"
                                    }
                                    = $seq_hash_ref {
                                        $fich1[
                                    0
                                    ]} {
                                        $variant_pos+$i
                                    } {
                                        "nt_pos"
                                    };
                                }
                                #}
								#}
                            }
                        }
                    }
                }
            }
            else {
                printf "##########ERROR MESSAGE : print $fich1[3]>$fich1[4] position $fich1[1] not simple mutation\n";
            }
        }
    }
    close IN;
    #print "$mut_nb1 mutation in hap1 and $mut_nb2 mutation in hap2 \n";
	#print Dumper (\%mut_hash);
	#print Dumper (\%mut_hash_pos);
    return (\%seq_hash_hap1, \%seq_hash_hap2, \%seq_hash_ref1, \%seq_hash_ref2, \%mut_hash, \%mut_hash_pos, $mut_nb1, $mut_nb2);
}
#tludwig is the following line normal  ?
1;
