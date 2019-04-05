package conf::config;
require Exporter;
use warnings;

@ISA = qw(Exporter);
@EXPORT = qw($genome_reference $CCDS_file $samtools $java $clin_var $vcftools );


# reference
############

$genome_reference="/PUBLIC_DATA/ReferenceGenomes/GRCh37.p5/reference.b37.fasta";

$CCDS_file="/PUBLIC_DATA/Annotation/CCDS/GRCh37/CCDS.current.txt";

$CCDS_NM_file="/PUBLIC_DATA/Annotation/CCDS/GRCh37/CCDS2Sequence.current.txt";

$clin_var="/PUBLIC_DATA/Annotation/Clinvar/clinvar_20180401.vcf.gz";

#bin list
#########

$samtools="/PROGS/EXTERN/samtools/samtools-1.3/samtools";
$vcftools="/PROGS/EXTERN/vcftools/vcftools_0.1.12b/bin/vcftools";
$java="java";
1;
